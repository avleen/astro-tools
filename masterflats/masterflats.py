#!/usr/bin/env python3
"""Turn each night's raw flats into one master flat per filter, next to the raw flats.

Walks the given folders (skipping any _selected* output trees), groups raw flats
by folder, night, filter, exposure, binning, gain and rotator angle, and writes

  <flat folder>/MasterFlat_<night>_<filter>_<exp>s.xisf

Each flat has a flat-dark subtracted, chosen in this order:
  1. flat-darks of the same night and exposure found in the tree (IMAGETYP DARK or DARKFLAT)
  2. a master dark of the same exposure from --darks (temperature within 2 C)
  3. a master bias from --darks
  4. the master dark nearest in exposure (dark current over a flat exposure is ~nil on
     a cooled CMOS; what has to come off is the camera offset)
It is then normalised to its median, and the flats are averaged with kappa-sigma
clipping (cosmic rays, satellites). The result is a Float32 XISF with
IMAGETYP = 'Master Flat', which WBPP uses as-is, zstd-compressed with byte shuffling
(as PixInsight's own "zstd + byte shuffling" save option; about 60% of the raw size).

Idempotent: the master records which flats and dark went into it. A master that is
up to date is left alone (one from an older version has its pixels compressed and its
header fixed in place, not rebuilt); new flats or a new flat-dark rebuild it. Once raw flats are
deleted (--delete-raw, or by hand) the master is kept. Safe to run every day.

subselect.py --copy picks the masters up automatically and copies them instead of
the raw flats.

Examples:
  masterflats.py "D:/Temp/Astro/Downloads" --darks "D:/Astro/Darks -10C/master" -n   # what would happen
  masterflats.py "D:/Temp/Astro/Downloads" --darks "D:/Astro/Darks -10C/master"      # build
  masterflats.py "D:/Temp/Astro/Downloads" --darks "D:/Astro/Darks -10C/master" --delete-raw

Requirements: Python 3.9+, numpy, zstandard (lz4 too, to read lz4-compressed flats).
"""
import argparse, html, json, os, re, struct, sys, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from xml.sax.saxutils import escape, quoteattr

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "subselect"))
from subselect import EXTS, _xisf_xml, expand, fix_argv, night_of, read_header, read_image, temp_of  # noqa: E402

VERSION = "1"
KEEP_KEYS = ("INSTRUME", "TELESCOP", "FOCALLEN", "XPIXSZ", "YPIXSZ", "XBINNING", "YBINNING", "GAIN", "OFFSET",
             "EGAIN", "SET-TEMP", "CCD-TEMP", "ROTATANG", "OBJECT", "DATE-OBS", "DATE-LOC")


def safe_name(s):
    return re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", s)


# ---------------------------------------------------------------- finding frames
def scan(paths):
    """Returns (raw flats, raw flat-darks) as dicts with path, file, hdr, night, filter, exp."""
    flats, fdarks = [], []
    for top in paths:
        if os.path.isfile(top):
            walk = [(os.path.dirname(top), [], [os.path.basename(top)])]
        elif os.path.isdir(top):
            walk = os.walk(top)
        else:
            print(f"  not found: {top}", file=sys.stderr)
            continue
        for root, dirs, files in walk:
            dirs[:] = sorted(d for d in dirs if not d.startswith("_selected"))    # subselect output: copies
            for fn in sorted(files):
                p = os.path.join(root, fn)
                if not fn.lower().endswith(EXTS):
                    continue
                try:
                    h = read_header(p)
                except Exception as e:
                    print(f"  skipping {p}: {e}", file=sys.stderr)
                    continue
                t = h.get("IMAGETYP", "").upper()
                if "MASTER" in t:
                    continue
                f = dict(path=p, file=fn, hdr=h, night=night_of(p, h), filter=h.get("FILTER", "") or "NOFILTER",
                         exp=float(h.get("EXPTIME", h.get("EXPOSURE", 0)) or 0))
                if "DARK" in t:
                    fdarks.append(f)
                elif "FLAT" in t:
                    flats.append(f)
    return flats, fdarks


def load_masters(ddirs):
    """Master darks and biases: list of (exposure s, temp C or None, path)."""
    out = []
    for ddir in ddirs:
        for root, _, files in os.walk(ddir):
            for fn in sorted(files):
                p = os.path.join(root, fn)
                if not fn.lower().endswith(EXTS):
                    continue
                try:
                    h = read_header(p)
                except Exception as e:
                    print(f"  skipping dark {p}: {e}", file=sys.stderr)
                    continue
                t = h.get("IMAGETYP", "").upper()
                if "MASTER" in t and ("DARK" in t or "BIAS" in t) and "FLAT" not in t:
                    out.append((float(h.get("EXPTIME", 0) or 0), temp_of(h, p), p))
    return out


def groups(flats):
    """Raw flats that belong in one master: same folder, night, filter, exposure, binning, gain, rotation."""
    g = defaultdict(list)
    for f in flats:
        h = f["hdr"]
        rot = h.get("ROTATANG", h.get("ROTATOR", ""))
        try:
            rot = f"{float(rot):.1f}"
        except ValueError:
            pass
        g[(os.path.dirname(f["path"]), f["night"], f["filter"], f["exp"], h.get("XBINNING", "1"),
           h.get("GAIN", ""), h.get("OFFSET", ""), rot)].append(f)
    out, names = [], defaultdict(int)
    for key in sorted(g, key=lambda k: min(f["path"] for f in g[k])):
        d, night, flt, exp = key[:4]
        name = safe_name(f"MasterFlat_{night}_{flt}_{exp:g}s")
        names[(d, name)] += 1
        if names[(d, name)] > 1:                     # two rotations/gains in one folder
            name += f"_{names[(d, name)]}"
        out.append(dict(night=night, filter=flt, exp=exp, gain=key[5], path=os.path.join(d, name + ".xisf"),
                        files=g[key]))
    return out


def pick_dark(grp, fdarks, masters):
    """(description, [paths]) of what to subtract from this group's flats, or None."""
    raw = [d["path"] for d in fdarks if d["night"] == grp["night"] and abs(d["exp"] - grp["exp"]) < 0.05
           and d["hdr"].get("GAIN", "") == grp["gain"]]
    if raw:
        return f"{len(raw)} flat-dark(s) from the same night", raw
    temp = temp_of(grp["files"][0]["hdr"])
    ok = [d for d in masters if temp is None or d[1] is None or abs(d[1] - temp) <= 2]
    exact = [d for d in ok if abs(d[0] - grp["exp"]) < 0.05]
    bias = [d for d in ok if d[0] < 0.01]
    if exact:
        return f"master dark {os.path.basename(exact[0][2])}", [exact[0][2]]
    if bias:
        return f"master bias {os.path.basename(bias[0][2])}", [bias[0][2]]
    if ok:
        d = min(ok, key=lambda d: abs(d[0] - grp["exp"]))
        return f"master dark {os.path.basename(d[2])} (nearest exposure)", [d[2]]
    return None


# ---------------------------------------------------------------- idempotency
def recorded_inputs(path):
    """The inputs recorded in a master we wrote, or None if there isn't one (or it isn't ours)."""
    try:
        with open(path, "rb") as f:
            xml = _xisf_xml(f)
    except (OSError, ValueError):
        return None
    m = re.search(r'<Property id="MasterFlats:Inputs" type="String">(.*?)</Property>', xml, re.S)
    return json.loads(html.unescape(m.group(1))) if m else None


def signature(grp, dark_paths):
    """By name and size, not path or mtime, so copies and moved trees still match."""
    return dict(version=VERSION, flats=sorted(f"{f['file']}|{os.path.getsize(f['path'])}" for f in grp["files"]),
                darks=sorted(f"{os.path.basename(p)}|{os.path.getsize(p)}" for p in dark_paths))


def up_to_date(prev, sig):
    if not prev or prev.get("version") != VERSION or not set(sig["flats"]) <= set(prev["flats"]):
        return False                                  # none, old, or there are flats it doesn't have
    all_present = set(prev["flats"]) == set(sig["flats"])
    return prev["darks"] == sig["darks"] or not all_present   # raw flats partly deleted: can't rebuild better


# ---------------------------------------------------------------- integration
def _as_uint16(img):
    u = img.astype(np.uint16)
    return u if np.array_equal(u, img) else img       # halves the memory for raw camera frames


def integrate(paths, dark, jobs, kappa=3.0, stripe=256):
    """Dark-subtract, normalise each flat to its median, then a kappa-sigma clipped mean per pixel.
    Returns the master in ADU, at the mean flat level."""
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        stack = list(ex.map(lambda p: _as_uint16(read_image(p)[0]), paths))
    if any(s.shape != stack[0].shape for s in stack) or dark.shape != stack[0].shape:
        raise ValueError("flats and flat-dark differ in size")
    scales = np.array([np.median(s[::4, ::4] - dark[::4, ::4]) for s in stack], np.float32)
    if (scales <= 0).any():
        raise ValueError("a flat is not brighter than its dark")
    out = np.empty(stack[0].shape, np.float32)
    for r0 in range(0, out.shape[0], stripe):
        x = np.stack([s[r0:r0 + stripe] for s in stack]).astype(np.float32)
        x -= dark[r0:r0 + stripe]
        x /= scales[:, None, None]
        med = np.median(x, 0)
        keep = np.abs(x - med) <= kappa * np.maximum(1.4826 * np.median(np.abs(x - med), 0), 1e-6)
        out[r0:r0 + stripe] = (x * keep).sum(0) / np.maximum(keep.sum(0), 1)
    return out * scales.mean()


CREATOR = f"astro-tools masterflats {VERSION}"

ZSTD_LEVEL = 9       # 59% of raw on a real master; level 19 gets 58% and takes 20x longer


def _compress(raw, item=4):
    """zstd + byte shuffle, as PixInsight writes it: returns (data, XISF compression attribute)."""
    import zstandard
    sh = np.frombuffer(raw, np.uint8).reshape(-1, item).T.tobytes()
    return (zstandard.ZstdCompressor(level=ZSTD_LEVEL, threads=-1).compress(sh),
            f"zstd+sh:{len(raw)}:{item}")


def _write_xisf(path, xml, data):
    """xml(offset) -> header bytes. Data block aligned to 4096 like PixInsight. Written atomically."""
    off = -(-(16 + len(xml(10 ** 12))) // 4096) * 4096
    x = xml(off)
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "wb") as f:
        f.write(b"XISF0100" + struct.pack("<II", len(x), 0) + x)
        f.write(b"\0" * (off - 16 - len(x)))
        f.write(data)
    os.replace(tmp, path)


def _metadata():
    """File-level properties: only the XISF: namespace is allowed here, and these two are mandatory
    (XISF 1.0 spec 11.4.1: a TimePoint, and the creating application's name and version)."""
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    return ('<Metadata>\n'
            f'<Property id="XISF:CreationTime" type="TimePoint" value="{now}"/>\n'
            f'<Property id="XISF:CreatorApplication" type="String">{CREATOR}</Property>\n'
            '</Metadata>\n')


INPUTS_RE = r'<Property id="MasterFlats:Inputs" type="String">.*?</Property>\n?'


def _stale(xml, tag):
    """What an older masterflats left out of date in this header: [] if nothing."""
    meta = re.search(r"<Metadata>.*?</Metadata>", xml, re.S)
    meta = meta.group(0) if meta else ""
    return ((["compressed"] if "compression=" not in tag else [])
            + (["header fixed"] if "MasterFlats:" in meta or 'XISF:CreationTime" type="TimePoint"' not in meta
               or f">{CREATOR}<" not in meta else []))


def refresh_in_place(path, dry_run=False):
    """Bring a master from an older masterflats up to date without rebuilding it: compress the pixels
    if they aren't, and move our inputs out of the file-level <Metadata> (PixInsight only allows XISF:
    properties there) and add XISF:CreationTime. Pixels are unchanged. Returns what was done, [] if nothing."""
    with open(path, "rb") as f:
        xml = _xisf_xml(f)
        tag = re.search(r"<Image\b[^>]*>", xml).group(0)
        todo = _stale(xml, tag)
        if dry_run or not todo:
            return todo
        off, size = map(int, re.search(r'location="attachment:(\d+):(\d+)"', tag).groups())
        f.seek(off)
        data = f.read(size)
    comp = re.search(r'compression="([^"]*)"', tag)
    if comp is None:
        data, comp = _compress(data)
    else:
        comp = comp.group(1)
    inputs = re.search(INPUTS_RE, xml, re.S).group(0).rstrip("\n") + "\n"
    body = re.sub(r"<Metadata>.*?</Metadata>\n?", "", re.sub(INPUTS_RE, "", xml, flags=re.S), flags=re.S)
    body = body.replace("</Image>", inputs + "</Image>", 1).replace("</xisf>", _metadata() + "</xisf>")
    plain = re.sub(r' compression="[^"]*"', "", tag)
    _write_xisf(path, lambda o: body.replace(tag, re.sub(
        r'location="attachment:\d+:\d+"', f'location="attachment:{o}:{len(data)}" compression="{comp}"', plain)
    ).encode(), data)
    return todo


def write_master(path, img, hdr, grp, sig, dark_desc, kappa):
    """Float32 [0,1] zstd+sh XISF that WBPP recognises as a master flat. Written atomically."""
    data, comp = _compress(np.ascontiguousarray(np.clip(img / 65535, 0, 1), "<f4").tobytes())
    h, w = img.shape

    def card(k, v, c=""):
        if not isinstance(v, (int, float)):
            try:
                float(v)
            except ValueError:
                v = "'" + str(v).replace("'", "''") + "'"
        return f"<FITSKeyword name={quoteattr(k)} value={quoteattr(str(v))} comment={quoteattr(c)}/>\n"
    cards = (card("IMAGETYP", "Master Flat", "Type of image") + card("FILTER", grp["filter"])
             + card("EXPTIME", grp["exp"], "[s] exposure of each flat") + card("EXPOSURE", grp["exp"])
             + card("NCOMBINE", len(grp["files"]), "number of flats combined")
             + "".join(card(k, hdr[k]) for k in KEEP_KEYS if hdr.get(k, "") != "")
             + card("HISTORY", "", f"masterflats: {len(grp['files'])} flats minus {dark_desc}, "
                                   f"median-normalised, {kappa:g}-sigma clipped mean"))
    inputs = f'<Property id="MasterFlats:Inputs" type="String">{escape(json.dumps(sig))}</Property>\n'
    xml = lambda off: ('<?xml version="1.0" encoding="UTF-8"?>\n'
                       '<xisf version="1.0" xmlns="http://www.pixinsight.com/xisf" '
                       'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                       'xsi:schemaLocation="http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd">\n'
                       f'<Image geometry="{w}:{h}:1" sampleFormat="Float32" bounds="0:1" imageType="MasterFlat" '
                       f'colorSpace="Gray" location="attachment:{off}:{len(data)}" compression="{comp}">\n'
                       f'{cards}{inputs}</Image>\n'
                       f'{_metadata()}</xisf>\n').encode()
    _write_xisf(path, xml, data)


# ---------------------------------------------------------------- main
def main(argv=None):
    for st in (sys.stdout, sys.stderr):              # cp1252 console / redirected file: never die on a character
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="folders (searched recursively) or flat files")
    ap.add_argument("--darks", action="append", default=[],
                    help="folder of master darks/biases, matched by exposure and temperature; repeatable")
    ap.add_argument("--delete-raw", action="store_true", help="delete raw flats once they are in a master")
    ap.add_argument("--min-flats", type=int, default=3, help="fewest flats to build a master from (default 3)")
    ap.add_argument("--kappa", type=float, default=3.0, help="sigma-clipping threshold (default 3)")
    ap.add_argument("-n", "--dry-run", action="store_true", help="show what would be built; change nothing")
    ap.add_argument("-j", "--jobs", type=int, default=min(4, os.cpu_count() or 1),
                    help="threads reading flats (default %(default)s)")
    a = ap.parse_args(fix_argv(sys.argv[1:] if argv is None else list(argv)))
    a.paths, a.darks = expand(a.paths), expand(a.darks)
    if not a.paths:
        sys.exit("no input folders")

    flats, fdarks = scan(a.paths)
    masters = load_masters(a.darks)
    print(f"{len(flats)} raw flats, {len(fdarks)} raw flat-darks, {len(masters)} master darks/biases")
    status = defaultdict(int)
    dark_cache = {}
    for grp in groups(flats):
        label = f"  NIGHT_{grp['night']}  {grp['filter']:>4} {grp['exp']:g}s {len(grp['files']):3d} flats"
        name = os.path.basename(grp["path"])
        dk = pick_dark(grp, fdarks, masters)
        if dk is None:
            print(f"{label}  NOT BUILT: no flat-dark, bias or master dark at this temperature (use --darks)")
            status["not built"] += 1
            continue
        sig = signature(grp, dk[1])
        prev = recorded_inputs(grp["path"])
        if prev is None and os.path.exists(grp["path"]):
            print(f"{label}  NOT BUILT: {name} exists and wasn't made by masterflats; leaving it alone")
            status["not built"] += 1
            continue
        if up_to_date(prev, sig):
            note = ""
            try:
                done = refresh_in_place(grp["path"], a.dry_run)
                note = done and f" ({'would be ' if a.dry_run else ''}{' and '.join(done)})"
            except Exception as e:
                note = f" (updating it FAILED: {type(e).__name__}: {e})"
                status["failed"] += 1
            print(f"{label}  up to date: {name}{note or ''}")
            status["up to date"] += 1
        elif len(grp["files"]) < a.min_flats:
            print(f"{label}  NOT BUILT: fewer than {a.min_flats} flats")
            status["not built"] += 1
            continue
        elif a.dry_run:
            print(f"{label}  would build {name} (minus {dk[0]})")
            status["to build"] += 1
            continue
        else:
            t0 = time.time()
            try:
                key = tuple(dk[1])
                if key not in dark_cache:
                    dark_cache.clear()                # one dark in memory at a time
                    ds = [read_image(p)[0] for p in key]
                    dark_cache[key] = np.median(np.stack(ds), 0) if len(ds) > 1 else ds[0]
                img = integrate([f["path"] for f in grp["files"]], dark_cache[key], a.jobs, a.kappa)
                write_master(grp["path"], img, grp["files"][0]["hdr"], grp, sig, dk[0], a.kappa)
            except Exception as e:
                print(f"{label}  FAILED: {type(e).__name__}: {e}")
                status["failed"] += 1
                continue
            print(f"{label}  built {name} (minus {dk[0]}, {time.time() - t0:.0f} s)")
            status["built"] += 1
        if a.delete_raw and not a.dry_run:
            done = set(recorded_inputs(grp["path"])["flats"])
            gone = [f for f in grp["files"] if f"{f['file']}|{os.path.getsize(f['path'])}" in done]
            for f in gone:
                os.remove(f["path"])
            print(f"{' ' * len(label)}  deleted {len(gone)} raw flats (all are in the master)")
    print(", ".join(f"{n} {k}" for k, n in status.items()) or "nothing to do")
    return 1 if status["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
