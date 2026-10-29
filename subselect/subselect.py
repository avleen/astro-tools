#!/usr/bin/env python3
"""Grade astrophotography light subs and copy the keepers into a stacking-ready tree.

Reads XISF (as written by N.I.N.A., including lz4/zstd + byte-shuffle compression)
and uncompressed FITS. For every light sub it measures:

  FWHM, minor axis   median star size, arcsec (Gaussian-weighted adaptive moments)
  ecc                median star eccentricity (0 = round; trailing/tilt raise it)
  trail              field-coherent elongation: a smear shared by all stars, i.e.
                     tracking/flexure, as opposed to random seeing/noise
  stars              detected star count (clouds, haze, bad focus drop it)
  bg, gradient       sky level and large-scale gradient (stray light, moon, cloud)
  sky e-, eff        optional (--darks): sky above dark in electrons and the
                     read-noise efficiency sky/(sky+RN^2)

Each sub is graded good / borderline / reject against the limits below. With
--copy, the keepers and the flats from the same night go to:

  OUT/NIGHT_<date>/LIGHT/...              good subs, all filters together
  OUT/NIGHT_<date>/FLAT/...               that night's flats, for the filters used (only the
                                          master flat, if there is one: see masterflats.py)
  OUT/NIGHT_<date>/DARKFLAT/...           that night's flat-darks, if any
  OUT_borderline/NIGHT_<date>/LIGHT/...   borderline subs

Add OUT as one folder in PixInsight WBPP (and OUT_borderline as well for a second,
comparison run). WBPP sorts files by type and filter itself. Add NIGHT as a
grouping keyword, Pre-processing only, so each night's lights use that night's
flats but all nights stack together.

Night = the evening date: taken from a NIGHT_YYYY-MM-DD folder in the source path
if there is one, else local time (DATE-LOC) minus 12 h.

Examples:
  subselect.py "D:/Astro/Heart Nebula"                      # report only
  subselect.py "Heart Nebula" --copy                        # copy to "Heart Nebula/_selected" (+ _selected_borderline)
  subselect.py lights/ --darks masters/ --csv report.csv    # + sky background, CSV
  subselect.py lights/ --fwhm-good 3.5 --ecc-good 0.6       # stricter limits

Requirements: Python 3.9+, numpy, scipy; lz4 and/or zstandard for compressed XISF.
"""
import argparse, csv, glob, hashlib, json, os, re, shutil, struct, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta

import numpy as np
from scipy import ndimage

np.seterr(all="ignore")

VERSION = "1"
CACHE = os.environ.get("SUBSELECT_CACHE") or (
    os.path.join(os.environ["LOCALAPPDATA"], "subselect") if os.name == "nt" and os.environ.get("LOCALAPPDATA")
    else os.path.join(os.path.expanduser("~"), ".cache", "subselect"))
SAT = 50000          # ADU: skip stars peaking above this (non-linear / clipped)
R = 12               # star cutout half-size, px
EXTS = (".xisf", ".fits", ".fit", ".fts")


# ---------------------------------------------------------------- file reading
def _xisf_xml(f):
    if f.read(8) != b"XISF0100":
        raise ValueError("not an XISF file")
    n = struct.unpack("<I", f.read(4))[0]
    f.read(4)
    return f.read(n).decode("utf-8", "replace")


def _xisf_keywords(xml):
    return {k: v.strip("'").strip() for k, v in
            re.findall(r'<FITSKeyword name="([^"]+)" value="([^"]*)"', xml)}


def _fits_header(f):
    hdr, size = {}, 0
    while True:
        block = f.read(2880)
        size += 2880
        if len(block) < 2880:
            raise ValueError("truncated FITS header")
        for i in range(0, 2880, 80):
            card = block[i:i + 80].decode("ascii", "replace")
            key = card[:8].strip()
            if key == "END":
                return hdr, size
            if card[8:10] == "= ":
                val = card[10:]
                if val.lstrip().startswith("'"):
                    m = re.match(r"\s*'((?:[^']|'')*)'", val)
                    val = m.group(1).replace("''", "'") if m else val
                else:
                    val = val.split("/")[0]
                hdr[key] = val.strip()


def read_header(path):
    with open(path, "rb") as f:
        if path.lower().endswith(".xisf"):
            return _xisf_keywords(_xisf_xml(f))
        return _fits_header(f)[0]


def read_image(path):
    """Returns (float32 image in ADU, header dict)."""
    with open(path, "rb") as f:
        if not path.lower().endswith(".xisf"):
            hdr, off = _fits_header(f)
            bp, w, h = int(hdr["BITPIX"]), int(hdr["NAXIS1"]), int(hdr["NAXIS2"])
            dt = {8: ">u1", 16: ">i2", 32: ">i4", -32: ">f4", -64: ">f8"}[bp]
            img = np.frombuffer(f.read(w * h * abs(bp) // 8), dt).reshape(h, w).astype(np.float32)
            img = img * float(hdr.get("BSCALE", 1)) + float(hdr.get("BZERO", 0))
            if bp < 0 and np.nanmax(img) <= 1.0:
                img *= 65535
            return img, hdr
        xml = _xisf_xml(f)
        tag = re.search(r"<Image\b[^>]*>", xml).group(0)
        w, h, _ = map(int, re.search(r'geometry="(\d+):(\d+):(\d+)"', tag).groups())
        off, size = map(int, re.search(r'location="attachment:(\d+):(\d+)"', tag).groups())
        fmt = re.search(r'sampleFormat="(\w+)"', tag).group(1)
        f.seek(off)
        data = f.read(size)
    comp = re.search(r'compression="([^"]+)"', tag)
    if comp:
        parts = comp.group(1).split(":")
        codec, usize = parts[0], int(parts[1])
        if codec.startswith("lz4"):
            import lz4.block
            data = lz4.block.decompress(data, uncompressed_size=usize)
        elif codec.startswith("zstd"):
            import zstandard
            data = zstandard.ZstdDecompressor().decompress(data, max_output_size=usize)
        elif codec.startswith("zlib"):
            import zlib
            data = zlib.decompress(data)
        else:
            raise ValueError(f"unsupported XISF compression {codec}")
        if codec.endswith("+sh"):
            isz = int(parts[2])
            data = np.frombuffer(data, np.uint8).reshape(isz, -1).T.tobytes()
    dt = {"UInt8": "u1", "UInt16": "<u2", "UInt32": "<u4", "Float32": "<f4", "Float64": "<f8"}[fmt]
    img = np.frombuffer(data, dt).reshape(h, w).astype(np.float32)
    if fmt.startswith("Float"):
        img *= 65535
    return img, _xisf_keywords(xml)


# ---------------------------------------------------------------- star measurement
def tile_medians(img, t=128):
    h, w = img.shape
    hh, ww = h // t * t, w // t * t
    tiles = img[:hh, :ww].reshape(hh // t, t, ww // t, t).transpose(0, 2, 1, 3).reshape(hh // t, ww // t, -1)
    return np.median(tiles, axis=2)


def adaptive_moments(cut, x0, y0):
    """Gaussian-weighted adaptive moments: centroid and covariance of the star."""
    yy, xx = np.mgrid[:cut.shape[0], :cut.shape[1]].astype(np.float64)
    cx, cy = x0, y0
    W = np.array([[4.0, 0.0], [0.0, 4.0]])
    for _ in range(30):
        try:
            Wi = np.linalg.inv(W)
        except np.linalg.LinAlgError:
            return None
        if not np.isfinite(Wi).all():
            return None
        dx, dy = xx - cx, yy - cy
        wt = cut * np.exp(-0.5 * (Wi[0, 0] * dx * dx + 2 * Wi[0, 1] * dx * dy + Wi[1, 1] * dy * dy))
        s = wt.sum()
        if s <= 0:
            return None
        ncx, ncy = (wt * xx).sum() / s, (wt * yy).sum() / s
        dx, dy = xx - ncx, yy - ncy
        M = np.array([[(wt * dx * dx).sum(), (wt * dx * dy).sum()],
                      [(wt * dx * dy).sum(), (wt * dy * dy).sum()]]) / s
        Wn = 2 * M
        if not np.isfinite(Wn).all() or Wn[0, 0] <= 0 or Wn[1, 1] <= 0 or np.linalg.det(Wn) <= 0 \
                or Wn[0, 0] > 60 or Wn[1, 1] > 60:
            return None
        conv = abs(ncx - cx) < 1e-3 and abs(ncy - cy) < 1e-3 and np.abs(Wn - W).max() < 1e-3
        cx, cy, W = ncx, ncy, Wn
        if conv:
            break
    return cx, cy, W


def find_stars(img):
    h, w = img.shape
    bgm = tile_medians(img)
    bg = ndimage.zoom(bgm, (h / bgm.shape[0], w / bgm.shape[1]), order=1)[:h, :w]
    d = img - bg
    noise = 1.4826 * np.median(np.abs(d[::4, ::4] - np.median(d[::4, ::4])))
    sm = ndimage.gaussian_filter(ndimage.median_filter(d, size=3), 1.5)   # median kills hot pixels
    sm_noise = 1.4826 * np.median(np.abs(sm[::4, ::4] - np.median(sm[::4, ::4])))
    peaks = (sm == ndimage.maximum_filter(sm, size=2 * R + 1)) & (sm > 10 * sm_noise)
    ys, xs = np.nonzero(peaks)
    keep = (xs > R + 2) & (xs < w - R - 3) & (ys > R + 2) & (ys < h - R - 3)
    stars = []
    for y, x in zip(ys[keep], xs[keep]):
        if img[y - R:y + R + 1, x - R:x + R + 1].max() > SAT:
            continue
        cut = d[y - R:y + R + 1, x - R:x + R + 1].astype(np.float64)
        m = adaptive_moments(cut, R, R)
        if m is None:
            continue
        cx, cy, C = m
        if not (np.isfinite(cx) and np.isfinite(cy)) or abs(cx - R) > 2 or abs(cy - R) > 2:
            continue
        ev = np.linalg.eigvalsh(C)
        if ev[0] < 0.6 ** 2:           # narrower than a pixel: hot pixel / cosmic ray
            continue
        flux = cut[cut > 0].sum()
        if flux / (noise * np.sqrt(np.pi * 4 * ev.mean()) + 1e-9) < 50:
            continue
        stars.append((x - R + cx, y - R + cy, C[0, 0], C[1, 1], C[0, 1]))
    return np.array(stars) if stars else np.zeros((0, 5)), bgm


def shape_stats(s):
    """Median star shape plus the field-coherent elongation (what a trail adds)."""
    if len(s) < 10:
        return dict(fwhm=np.nan, minor=np.nan, ecc=np.nan, trail=np.nan, pa=np.nan)
    Cxx, Cyy, Cxy = s[:, 2], s[:, 3], s[:, 4]
    tr, det = Cxx + Cyy, Cxx * Cyy - Cxy ** 2
    disc = np.sqrt(np.maximum((tr / 2) ** 2 - det, 0))
    a2, b2 = tr / 2 + disc, tr / 2 - disc
    # Seeing/noise elongations point every which way and cancel in the median
    # shear; a common trail of length L adds L^2/12 along its direction.
    d1, d2 = np.median(Cxx - Cyy), np.median(2 * Cxy)
    return dict(fwhm=2.3548 * np.median(np.sqrt(np.sqrt(a2 * b2))),
                minor=2.3548 * np.median(np.sqrt(b2)),
                ecc=float(np.median(np.sqrt(1 - b2 / a2))),
                trail=float(np.sqrt(12 * np.hypot(d1, d2))),
                pa=float(np.degrees(0.5 * np.arctan2(d2, d1)) % 180))


def measure(path):
    """Per-sub metrics in pixels/ADU, cached by path+size+mtime."""
    st = os.stat(path)
    key = hashlib.sha1(f"{os.path.normcase(os.path.abspath(path))}|{st.st_size}|{st.st_mtime}|{VERSION}".encode()).hexdigest()[:20]
    cpath = os.path.join(CACHE, key + ".json")
    try:
        with open(cpath, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        pass
    img, hdr = read_image(path)
    stars, bgm = find_stars(img)
    g = np.array([[np.median(b) for b in np.array_split(row, 3, axis=1)]
                  for row in np.array_split(bgm, 3, axis=0)])
    m = dict(shape_stats(stars), n=len(stars),
             bg=float(np.median(bgm)),
             bg_dark=float(np.median(bgm[bgm <= np.percentile(bgm, 10)])),   # sky, avoiding nebula
             gradient=float(g.max() - g.min()), hdr=hdr)
    m = {k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in m.items()}
    os.makedirs(CACHE, exist_ok=True)
    tmp = f"{cpath}.{os.getpid()}.tmp"            # atomic: parallel runs never see half a file
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(m, f)
    os.replace(tmp, cpath)
    return m


def _measure_safe(path):
    try:
        return path, measure(path), None
    except Exception as e:                      # one bad file shouldn't stop the run
        return path, None, f"{type(e).__name__}: {e}"


# ---------------------------------------------------------------- metadata
def night_of(path, hdr):
    m = re.search(r"NIGHT_(\d{4}-\d{2}-\d{2})", path.replace("\\", "/"))
    if m:
        return m.group(1)
    for k in ("DATE-LOC", "DATE-OBS"):
        if hdr.get(k):
            return (datetime.fromisoformat(hdr[k][:19]) - timedelta(hours=12)).date().isoformat()
    return "unknown"


def time_of(hdr):
    v = hdr.get("DATE-LOC") or hdr.get("DATE-OBS")
    return datetime.fromisoformat(v[:19]) if v else None


def pixel_scale(hdr, override):
    if override:
        return override
    try:
        return 206.265 * float(hdr["XPIXSZ"]) / float(hdr["FOCALLEN"])   # XPIXSZ already includes binning
    except (KeyError, ValueError, ZeroDivisionError):
        return None


def temp_of(hdr, path=""):
    for k in ("SET-TEMP", "CCD-TEMP"):
        try:
            return float(hdr[k])
        except (KeyError, ValueError):
            pass
    m = re.search(r"(?<![\d.])(-?\d+(?:\.\d+)?) ?C\b", path.replace("\\", "/").rpartition("/")[0])  # e.g. "Darks - -10C, Gain 100"
    return float(m.group(1)) if m else None


def load_darks(ddirs):
    """Master darks: list of (exposure s, sensor temp C or None, median ADU)."""
    out = []
    for ddir in ddirs:
        for p in sorted(glob.glob(os.path.join(glob.escape(ddir), "**", "*"), recursive=True)):
            if not p.lower().endswith(EXTS):
                continue
            try:
                h = read_header(p)
                t = h.get("IMAGETYP", "").upper()
                if "DARK" not in t or "FLAT" in t:
                    continue
                img, _ = read_image(p)
                out.append((float(h["EXPTIME"]), temp_of(h, p), float(np.median(img[::4, ::4]))))
            except Exception as e:
                print(f"  skipping dark {p}: {e}", file=sys.stderr)
    return out


def match_dark(darks, exp, temp):
    best = None
    for e, t, med in darks:
        if abs(e - exp) > 0.05:
            continue
        if temp is not None and t is not None and abs(t - temp) > 2:
            continue
        if best is None or (t is not None and temp is not None and abs(t - temp) < abs(best[1] - temp)):
            best = (med, t)
    return best[0] if best else None


# ---------------------------------------------------------------- grading
def grade(s, a):
    """Returns (tier, reasons). tier: good / borderline / reject."""
    bad, meh = [], []
    def check(val, good_lim, max_lim, label, fmt, higher_is_worse=True):
        if val is None:
            return
        worse = (lambda v, l: v > l) if higher_is_worse else (lambda v, l: v < l)
        if worse(val, max_lim):
            bad.append(f"{label} {val:{fmt}}")
        elif worse(val, good_lim):
            meh.append(f"{label} {val:{fmt}}")

    if s["n"] < a.min_stars or s["fwhm"] is None:
        return "reject", [f"only {s['n']} stars"]
    check(s["fwhm_as"], a.fwhm_good, a.fwhm_max, "FWHM", ".2f")
    check(s["ecc"], a.ecc_good, a.ecc_max, "ecc", ".2f")
    check(s["star_ratio"], a.stars_good, a.stars_min, "stars x", ".2f", higher_is_worse=False)
    check(s["bg_ratio"], a.bg_good, a.bg_max, "background x", ".2f")
    check(s["grad_ratio"], a.grad_good, a.grad_max, "gradient x", ".1f")
    return ("reject" if bad else "borderline" if meh else "good"), bad + meh


def median_by(subs, keyf, val):
    groups = defaultdict(list)
    for s in subs:
        if s[val] is not None:
            groups[keyf(s)].append(s[val])
    return {k: float(np.median(v)) for k, v in groups.items()}


# ---------------------------------------------------------------- main
def fix_argv(argv):
    """cmd.exe/PowerShell: "D:\\dir\\" makes \\" a literal quote that swallows the following args."""
    if os.name != "nt":
        return argv
    out = []
    for x in argv:
        head, q, rest = x.partition('"')
        out += [head] + rest.split() if q else [x]
    return [x for x in out if x]


def expand(paths):
    """Expand wildcards ourselves: cmd.exe and PowerShell pass them through unexpanded."""
    out = []
    for p in paths:
        if not glob.has_magic(p) or os.path.exists(p):
            out.append(p)
            continue
        hits = sorted(glob.glob(p, recursive=True))
        if not hits:
            print(f"  no match for {p}", file=sys.stderr)
        out += hits
    return out


def find_files(paths, exclude):
    norm = lambda p: os.path.normcase(os.path.abspath(p))
    exclude = {norm(p) for p in exclude}
    out = []
    for p in paths:
        if os.path.isfile(p):
            out.append(os.path.abspath(p))
            continue
        if not os.path.isdir(p):
            print(f"  not found: {p}", file=sys.stderr)
            continue
        for root, dirs, files in os.walk(p):
            dirs[:] = [d for d in dirs if norm(os.path.join(root, d)) not in exclude
                       and not d.startswith("_selected")]      # other targets' output trees
            out += [os.path.abspath(os.path.join(root, f)) for f in files if f.lower().endswith(EXTS)]
    return sorted(set(out), key=lambda p: (os.path.normcase(p), p))


def main(argv=None):
    for st in (sys.stdout, sys.stderr):           # cp1252 console / redirected file: never die on a character
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="light files or folders (searched recursively)")
    ap.add_argument("--out", help="destination root (default: <first folder>/_selected); "
                                  "borderline subs go to <out>_borderline")
    ap.add_argument("--copy", action="store_true", help="copy good + borderline subs (default: report only)")
    ap.add_argument("--good-only", action="store_true", help="with --copy, skip borderline")
    ap.add_argument("--csv", help="write the per-sub table to this CSV file")
    ap.add_argument("--darks", action="append", default=[],
                    help="folder of master darks (matched by exposure and temperature) to report sky e- and "
                         "read-noise efficiency; can be given more than once")
    ap.add_argument("--read-noise", type=float, default=1.4, help="camera read noise, e- (default 1.4)")
    ap.add_argument("--gain", type=float, help="e-/ADU if the header has no EGAIN")
    ap.add_argument("--scale", type=float, help="arcsec/px if the header has no XPIXSZ/FOCALLEN")
    ap.add_argument("-j", "--jobs", type=int, default=min(4, os.cpu_count() or 1),
                    help="parallel workers (each needs ~1 GB RAM for a 26 MP frame; default %(default)s)")
    g = ap.add_argument_group("limits (good / reject-beyond; in between = borderline)")
    g.add_argument("--fwhm-good", type=float, default=4.5, help="arcsec (default %(default)s)")
    g.add_argument("--fwhm-max", type=float, default=6.0, help="arcsec (default %(default)s)")
    g.add_argument("--ecc-good", type=float, default=0.72, help="(default %(default)s)")
    g.add_argument("--ecc-max", type=float, default=0.82, help="(default %(default)s)")
    g.add_argument("--stars-good", type=float, default=0.5,
                   help="star count vs median for that filter+exposure (default %(default)s)")
    g.add_argument("--stars-min", type=float, default=0.25, help="(default %(default)s)")
    g.add_argument("--min-stars", type=int, default=30, help="absolute minimum star count (default %(default)s)")
    g.add_argument("--bg-good", type=float, default=1.5,
                   help="background vs median for that filter+exposure+night (default %(default)s)")
    g.add_argument("--bg-max", type=float, default=2.0, help="(default %(default)s)")
    g.add_argument("--grad-good", type=float, default=2.0,
                   help="gradient vs median for that filter+exposure+night (default %(default)s)")
    g.add_argument("--grad-max", type=float, default=3.0, help="(default %(default)s)")
    a = ap.parse_args(fix_argv(sys.argv[1:] if argv is None else list(argv)))
    a.paths, a.darks = expand(a.paths), expand(a.darks)
    if not a.paths:
        sys.exit("no input files or folders")
    if os.name == "nt":
        a.jobs = min(a.jobs, 61)                  # ProcessPoolExecutor limit on Windows
    a.jobs = max(1, a.jobs)

    first_dir = next((p for p in a.paths if os.path.isdir(p)), os.path.dirname(os.path.abspath(a.paths[0])))
    out_root = os.path.abspath(a.out or os.path.join(first_dir, "_selected"))   # abspath drops a trailing slash
    out_border = out_root + "_borderline"

    files, flats, masters = [], [], defaultdict(list)
    for p in find_files(a.paths, {out_root, out_border}):
        try:
            h = read_header(p)
        except Exception as e:
            print(f"  skipping {p}: {e}", file=sys.stderr)
            continue
        t = h.get("IMAGETYP", "LIGHT").upper()
        if "MASTER" in t:
            if "FLAT" in t and "DARK" not in t:
                masters[(night_of(p, h), h.get("FILTER", "") or "NOFILTER")].append(p)
        elif "FLAT" in t:
            flats.append(dict(path=p, file=os.path.basename(p), filter=h.get("FILTER", "") or "NOFILTER",
                              night=night_of(p, h), kind="DARKFLAT" if "DARK" in t else "FLAT"))
        elif "LIGHT" in t:
            files.append(p)
    if not files:
        sys.exit("no light frames found")
    print(f"measuring {len(files)} light subs with {a.jobs} worker(s) (cached results are instant)...", file=sys.stderr)

    subs = []
    with ProcessPoolExecutor(max_workers=min(a.jobs, len(files))) as ex:
        for i, (p, m, err) in enumerate(ex.map(_measure_safe, files, chunksize=1), 1):
            print(f"\r  {i}/{len(files)}", end="", file=sys.stderr, flush=True)
            if err:
                print(f"\n  FAILED {p}: {err}", file=sys.stderr)
                continue
            h = m.pop("hdr")
            sc = pixel_scale(h, a.scale)
            s = dict(m, path=p, file=os.path.basename(p), hdr=h,
                     filter=h.get("FILTER", "") or "NOFILTER",
                     exp=float(h.get("EXPTIME", h.get("EXPOSURE", 0)) or 0),
                     night=night_of(p, h), time=time_of(h), scale=sc)
            for k in ("fwhm", "minor", "trail"):
                s[k + "_as"] = s[k] * sc if (s[k] is not None and sc) else None
            subs.append(s)
    print(file=sys.stderr)

    # relative checks: compare against similar subs
    fe = lambda s: (s["filter"], s["exp"])
    fen = lambda s: (s["filter"], s["exp"], s["night"])
    med_n, med_bg, med_gr = median_by(subs, fe, "n"), median_by(subs, fen, "bg"), median_by(subs, fen, "gradient")
    for s in subs:
        s["star_ratio"] = s["n"] / med_n[fe(s)] if med_n.get(fe(s)) else None
        s["bg_ratio"] = s["bg"] / med_bg[fen(s)] if med_bg.get(fen(s)) else None
        s["grad_ratio"] = s["gradient"] / med_gr[fen(s)] if med_gr.get(fen(s)) else None
        if s["fwhm_as"] is None and s["fwhm"] is not None:
            s["fwhm_as"] = s["fwhm"]                     # no scale known: limits apply in px
        s["tier"], s["why"] = grade(s, a)

    darks = load_darks(a.darks) if a.darks else {}
    for s in subs:
        s["sky_e"] = s["eff"] = None
        if darks and s["exp"]:
            d = match_dark(darks, s["exp"], temp_of(s["hdr"]))
            gain = float(s["hdr"].get("EGAIN", 0) or 0) or a.gain
            if d is not None and gain:
                s["sky_e"] = (s["bg_dark"] - d) * gain
                s["eff"] = s["sky_e"] / (s["sky_e"] + a.read_noise ** 2) * 100

    subs.sort(key=lambda s: (s["night"], s["time"] or datetime.min, s["file"]))

    # ---- table
    print(f"\n{'night':10} {'local':8} {'flt':>4} {'exp':>5} {'stars':>6} {'FWHM':>5} {'minor':>5} "
          f"{'ecc':>5} {'trail':>5} {'bg x':>5} {'grad x':>6}" + ("  sky e-  eff" if darks else "") + "  grade")
    f2 = lambda v, w=5, p=2: f"{v:{w}.{p}f}" if v is not None else " " * (w - 1) + "-"
    for s in subs:
        line = (f"{s['night']:10} {s['time'].strftime('%H:%M:%S') if s['time'] else '?':8} {s['filter']:>4} "
                f"{s['exp']:5.0f} {s['n']:6d} {f2(s['fwhm_as'])} {f2(s['minor_as'])} {f2(s['ecc'])} "
                f"{f2(s['trail_as'])} {f2(s['bg_ratio'])} {f2(s['grad_ratio'], 6, 1)}")
        if darks:
            line += f"  {f2(s['sky_e'], 6, 1)} {f2(s['eff'], 4, 0)}%"
        tag = {"good": "good", "borderline": "BORDERLINE", "reject": "REJECT"}[s["tier"]]
        print(line + f"  {tag}" + (f"  ({', '.join(s['why'])})" if s["why"] else ""))
    sc = next((s["scale"] for s in subs if s["scale"]), None)
    print(f"\nFWHM/minor/trail in arcsec ({sc:.3f}\"/px)" if sc else "\nFWHM/minor/trail in PIXELS (no scale: use --scale)")
    print("bg x / grad x = relative to the median of the same filter, exposure and night; "
          "stars compared within filter + exposure")

    # ---- totals
    tot = defaultdict(lambda: defaultdict(float))
    cnt = defaultdict(lambda: defaultdict(int))
    for s in subs:
        tot[s["filter"]][s["tier"]] += s["exp"] / 60
        cnt[s["filter"]][s["tier"]] += 1
    print(f"\n{'filter':>6}  {'good':>14}  {'borderline':>14}  {'reject':>14}")
    for flt in sorted(tot):
        cells = [f"{cnt[flt][t]:3d} = {tot[flt][t]:6.1f} min" for t in ("good", "borderline", "reject")]
        print(f"{flt:>6}  " + "  ".join(cells))

    if a.csv:
        cols = ["night", "time", "filter", "exp", "file", "n", "fwhm_as", "minor_as", "ecc", "trail_as", "pa",
                "bg", "bg_ratio", "gradient", "grad_ratio", "star_ratio", "sky_e", "eff", "tier", "why", "path"]
        with open(a.csv, "w", newline="", encoding="utf-8-sig") as f:      # -sig: Excel reads it as UTF-8
            w = csv.writer(f)
            w.writerow(cols)
            for s in subs:
                w.writerow([("; ".join(s[c]) if c == "why" else s[c]) for c in cols])
        print(f"\nwrote {a.csv}")

    if not a.copy:
        print(f"\nreport only. Add --copy to copy good + borderline subs to {out_root}")
        return
    tiers = ("good",) if a.good_only else ("good", "borderline")
    keep = [s for s in subs if s["tier"] in tiers]
    used = {(s["night"], s["filter"]) for s in keep}
    flat_sets = {(f["night"], f["filter"]) for f in flats if f["kind"] == "FLAT"} | set(masters)
    plan = [(s["path"], os.path.join(out_root if s["tier"] == "good" else out_border,
                                     f"NIGHT_{s['night']}", "LIGHT", s["file"])) for s in keep]
    # a master flat (masterflats.py) stands in for that night+filter's raw flats and flat-darks
    plan += [(m, os.path.join(out_root, f"NIGHT_{k[0]}", "FLAT", os.path.basename(m)))
             for k in sorted(used & set(masters)) for m in masters[k]]
    raw_used = used - set(masters)
    flat_nights = {n for n, _ in raw_used}
    plan += [(f["path"], os.path.join(out_root, f"NIGHT_{f['night']}", f["kind"], f["file"])) for f in flats
             if (f["night"], f["filter"]) in raw_used or (f["kind"] == "DARKFLAT" and f["night"] in flat_nights)]
    n = 0
    for src, dst in plan:
        if os.path.isfile(dst):
            ss, ds = os.stat(src), os.stat(dst)
            if ss.st_size == ds.st_size and abs(ss.st_mtime - ds.st_mtime) < 2:   # 2 s: FAT/exFAT mtime resolution
                continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        n += 1
    print(f"\ncopied {n} new file(s) ({len(plan)} planned) into {out_root}"
          + ("" if a.good_only else f" and {out_border}") + " - originals untouched")
    stale = defaultdict(int)                          # raw flats an earlier run copied, now beside a master
    for night, flt in used & set(masters):
        d = os.path.join(out_root, f"NIGHT_{night}", "FLAT")
        for p in glob.glob(os.path.join(glob.escape(d), "*")):
            h = read_header(p) if p.lower().endswith(EXTS) else {}
            if h and "MASTER" not in h.get("IMAGETYP", "").upper() and h.get("FILTER", "") == flt:
                stale[(night, flt)] += 1
    if stale:
        print("\nNOTE: these FLAT folders hold raw flats from an earlier run as well as the master flat. "
              "Delete the raw ones so WBPP uses only the master:")
        for (night, flt), c in sorted(stale.items()):
            print(f"  NIGHT_{night}/FLAT  {flt}: {c} raw flats")
    missing = sorted(used - flat_sets)
    if missing:
        print("\nWARNING: no flats found for these nights/filters (WBPP will need flats from another night):")
        for night, flt in missing:
            print(f"  NIGHT_{night}  {flt}")


if __name__ == "__main__":
    main()
