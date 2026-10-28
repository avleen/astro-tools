#!/usr/bin/env python3
"""Measure differential flexure from a set of light subs (XISF, as written by N.I.N.A.).

Two independent measurements:
  1. Within each sub: star trailing, from the field-coherent part of the star
     second moments. Seeing blurs stars round; a trail adds L^2/12 along its axis.
  2. Between subs: image drift (shift + rotation) from matched stars.
     Dithering must be off for this to mean anything.

If both agree (drift between subs predicts the trail within a sub, same
direction), the guider is holding the guide star while the main image moves:
that's differential flexure.

Usage:
  flextest.py FILE...                       # one sequence
  flextest.py --off F1 [F2...] --on F3 F4   # compare compensator off vs on
"""
import argparse, glob, hashlib, html, json, os, re, struct, sys
from datetime import datetime
import numpy as np
np.seterr(all="ignore")
from scipy import ndimage
from scipy.spatial import cKDTree

CACHE = os.environ.get("FLEXTEST_CACHE") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
SAT = 50000          # ADU; skip stars peaking above this (non-linear / clipped)
R = 12               # cutout half-size, px


# ---------------------------------------------------------------- console / args
def fix_console():
    """Don't die on characters the console or a redirect can't encode (cp1252 on Windows)."""
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


def expand(args):
    """Expand wildcards ourselves: cmd.exe and PowerShell pass them through literally."""
    out = []
    for a in args:
        if os.path.exists(a) or not any(c in a for c in "*?["):
            out.append(a); continue
        m = sorted(glob.glob(a, recursive=True))
        if not m:
            sys.exit(f"error: no files match {a}")
        out += m
    return out


# ---------------------------------------------------------------- XISF reading
def read_xisf(path):
    with open(path, "rb") as f:
        if f.read(8) != b"XISF0100":
            raise ValueError(f"{path}: not an XISF file")
        n = struct.unpack("<I", f.read(4))[0]; f.read(4)
        xml = f.read(n).decode("utf-8", "replace")
        tag = re.search(r"<Image\b[^>]*>", xml).group(0)
        w, h, _ = map(int, re.search(r'geometry="(\d+):(\d+):(\d+)"', tag).groups())
        off, size = map(int, re.search(r'location="attachment:(\d+):(\d+)"', tag).groups())
        fmt = re.search(r'sampleFormat="(\w+)"', tag).group(1)
        f.seek(off); data = f.read(size)
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
            raise ValueError(f"unsupported compression {codec}")
        if codec.endswith("+sh"):
            isz = int(parts[2])
            data = np.frombuffer(data, np.uint8).reshape(isz, -1).T.tobytes()
    dt = {"UInt16": "<u2", "Float32": "<f4", "UInt32": "<u4"}[fmt]
    img = np.frombuffer(data, dt).reshape(h, w).astype(np.float32)
    if fmt == "Float32":
        img *= 65535
    hdr = {k: html.unescape(v).strip("'").strip() for k, v in
           re.findall(r'<FITSKeyword name="([^"]+)" value="([^"]*)"', xml)}
    return img, hdr


# ---------------------------------------------------------------- star measurement
def background(img, t=128):
    h, w = img.shape
    hh, ww = h // t * t, w // t * t
    tiles = img[:hh, :ww].reshape(hh // t, t, ww // t, t).transpose(0, 2, 1, 3).reshape(hh // t, ww // t, -1)
    med = np.median(tiles, axis=2)
    return ndimage.zoom(med, (h / med.shape[0], w / med.shape[1]), order=1)[:h, :w]


def adaptive_moments(cut, x0, y0):
    """Gaussian-weighted adaptive moments (weight matched to the star). Returns
    centroid and covariance of the underlying (Gaussian-equivalent) star."""
    yy, xx = np.mgrid[:cut.shape[0], :cut.shape[1]].astype(np.float64)
    cx, cy = x0, y0
    W = np.array([[4.0, 0.0], [0.0, 4.0]])
    for _ in range(30):
        try:
            Wi = np.linalg.inv(W)
            if not np.isfinite(Wi).all():
                return None
        except np.linalg.LinAlgError:
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


_cache_warned = False


def cache_key(path):
    """From the file header (which holds the data checksum) and size, not the path or
    mtime, so the same file hits the same entry from Windows, WSL or a copy."""
    with open(path, "rb") as f:
        head = f.read(65536)
    return hashlib.sha1(head + f"|{os.path.getsize(path)}|v3".encode()).hexdigest()[:16]


def measure(path):
    global _cache_warned
    cpath = os.path.join(CACHE, cache_key(path) + ".npz")
    if os.path.exists(cpath):
        try:
            with np.load(cpath) as z:
                return z["stars"], json.loads(str(z["hdr"])), tuple(int(v) for v in z["shape"])
        except Exception:
            pass                       # corrupt / old format: re-measure

    img, hdr = read_xisf(path)
    h, w = img.shape
    d = img - background(img)
    noise = 1.4826 * np.median(np.abs(d[::4, ::4] - np.median(d[::4, ::4])))
    sm = ndimage.gaussian_filter(ndimage.median_filter(d, size=3), 1.5)   # median kills hot pixels
    sm_noise = 1.4826 * np.median(np.abs(sm[::4, ::4] - np.median(sm[::4, ::4])))
    peaks = (sm == ndimage.maximum_filter(sm, size=2 * R + 1)) & (sm > 10 * sm_noise)
    ys, xs = np.nonzero(peaks)
    keep = (xs > R + 2) & (xs < w - R - 3) & (ys > R + 2) & (ys < h - R - 3)
    ys, xs = ys[keep], xs[keep]

    stars = []
    for y, x in zip(ys, xs):
        cut = d[y - R:y + R + 1, x - R:x + R + 1].astype(np.float64)
        if img[y - R:y + R + 1, x - R:x + R + 1].max() > SAT:
            continue
        m = adaptive_moments(cut, R, R)
        if m is None:
            continue
        cx, cy, C = m
        if abs(cx - R) > 2 or abs(cy - R) > 2:
            continue
        ev = np.linalg.eigvalsh(C)
        if ev[0] < 0.6 ** 2:           # narrower than a pixel: hot pixel / cosmic ray
            continue
        if not (np.isfinite(cx) and np.isfinite(cy)):
            continue
        flux = cut[cut > 0].sum()
        snr = flux / (noise * np.sqrt(np.pi * 4 * ev.mean()) + 1e-9)
        if snr < 50:
            continue
        stars.append((x - R + cx, y - R + cy, C[0, 0], C[1, 1], C[0, 1], flux, snr))
    stars = np.array(stars) if stars else np.zeros((0, 7))

    try:
        os.makedirs(CACHE, exist_ok=True)
        tmp = f"{cpath}.{os.getpid()}.tmp"
        with open(tmp, "wb") as fh:
            np.savez(fh, stars=stars, hdr=np.array(json.dumps(hdr)), shape=np.array([h, w]))
        os.replace(tmp, cpath)
    except OSError as e:
        if not _cache_warned:
            print(f"warning: can't write cache in {CACHE} ({e}); set FLEXTEST_CACHE to a writable dir",
                  file=sys.stderr)
            _cache_warned = True
    return stars, hdr, (h, w)


# ---------------------------------------------------------------- per-frame metrics
def shape_stats(s):
    """Median shape, plus the field-coherent elongation (what a trail adds)."""
    Cxx, Cyy, Cxy = s[:, 2], s[:, 3], s[:, 4]
    tr = Cxx + Cyy
    det = Cxx * Cyy - Cxy ** 2
    disc = np.sqrt(np.maximum((tr / 2) ** 2 - det, 0))
    a2, b2 = tr / 2 + disc, tr / 2 - disc
    ecc = np.sqrt(1 - b2 / a2)
    # Coherent part: median shear components. Random (seeing/noise) cancels; a
    # common trail of length L adds (L^2/12) along its direction.
    d1, d2 = np.median(Cxx - Cyy), np.median(2 * Cxy)
    mag = np.hypot(d1, d2)
    return dict(fwhm=2.3548 * np.median(np.sqrt(np.sqrt(a2 * b2))),
                fwhm_min=2.3548 * np.median(np.sqrt(b2)),
                ecc=np.median(ecc),
                trail=np.sqrt(12 * mag),
                pa=np.degrees(0.5 * np.arctan2(d2, d1)) % 180,
                n=len(s))


def field_map(s, shape, n=3):
    h, w = shape
    out = []
    for gy in range(n):
        row = []
        for gx in range(n):
            m = ((s[:, 0] * n // w) == gx) & ((s[:, 1] * n // h) == gy)
            row.append(shape_stats(s[m]) if m.sum() >= 20 else None)
        out.append(row)
    return out


# ---------------------------------------------------------------- registration
def register(a, b, shape, max_shift=400):
    """Similarity transform mapping star positions in a -> b, about image centre."""
    h, w = shape
    ca, cb = np.array([w / 2, h / 2]), np.array([w / 2, h / 2])
    A = a[np.argsort(-a[:, 5])][:150, :2]
    B = b[np.argsort(-b[:, 5])][:150, :2]
    dxy = (B[None, :, :] - A[:, None, :]).reshape(-1, 2)
    dxy = dxy[(np.abs(dxy) < max_shift).all(1)]
    if len(dxy) == 0:
        return None
    H, ex, ey = np.histogram2d(dxy[:, 0], dxy[:, 1], bins=int(2 * max_shift / 2),
                               range=[[-max_shift, max_shift]] * 2)
    i, j = np.unravel_index(H.argmax(), H.shape)
    t = np.array([(ex[i] + ex[i + 1]) / 2, (ey[j] + ey[j + 1]) / 2])
    pa, pb = a[:, :2] - ca, b[:, :2] - cb
    tree = cKDTree(pb)
    p = np.array([1.0, 0.0, t[0], t[1]])     # x' = p0 x - p1 y + p2 ; y' = p1 x + p0 y + p3
    tol = 4.0
    for _ in range(6):
        pred = np.c_[p[0] * pa[:, 0] - p[1] * pa[:, 1] + p[2], p[1] * pa[:, 0] + p[0] * pa[:, 1] + p[3]]
        dist, idx = tree.query(pred, distance_upper_bound=tol)
        ok = np.isfinite(dist)
        if ok.sum() < 10:
            return None
        x, y = pa[ok, 0], pa[ok, 1]
        u, v = pb[idx[ok], 0], pb[idx[ok], 1]
        M = np.r_[np.c_[x, -y, np.ones_like(x), np.zeros_like(x)],
                  np.c_[y, x, np.zeros_like(x), np.ones_like(x)]]
        p = np.linalg.lstsq(M, np.r_[u, v], rcond=None)[0]
        res = np.r_[u, v] - M @ p
        tol = max(1.0, 4 * 1.4826 * np.median(np.abs(res)))
    return dict(dx=p[2], dy=p[3], rot=np.degrees(np.arctan2(p[1], p[0])) * 60,
                scale_ppm=(np.hypot(p[0], p[1]) - 1) * 1e6,
                rms=np.sqrt(np.mean(res ** 2)), nmatch=int(ok.sum()))


# ---------------------------------------------------------------- reporting
def t_of(hdr, key="DATE-OBS"):
    v = hdr.get(key) or hdr.get("DATE-OBS")
    if not v:
        raise ValueError("no DATE-OBS in header")
    return datetime.fromisoformat(v[:26].rstrip("Z"))


def warn(fr, prev):
    w = []
    if fr["st"]["trail"] > R * 1.5:
        w.append("trail longer than star cutout: value is a LOWER BOUND")
    if fr["st"]["n"] < 200:
        w.append("FEW STARS (cloud/trailed?)")
    if prev is not None:
        dra = (float(fr["hdr"].get("RA", 0)) - float(prev["hdr"].get("RA", 0))) * np.cos(np.radians(float(fr["hdr"].get("DEC", 0))))
        dt = max((fr["t"] - prev["t"]).total_seconds() / 60, 1e-3)
        if abs(dra) > 0.05 and abs(dra) / dt > 0.1:      # sidereal is 0.25 deg/min x cos(dec)
            w.append(f"MOUNT RA moved {dra:+.2f} deg in {dt:.0f} min - tracking stopped?")
        if prev["hdr"].get("FOCPOS") == fr["hdr"].get("FOCPOS") and \
                abs(float(fr["hdr"].get("FOCTEMP", 0)) - float(prev.get("focus_t0", prev["hdr"].get("FOCTEMP", 0)))) >= 2:
            w.append("focus unchanged over >=2C")
        else:
            fr["focus_t0"] = prev.get("focus_t0", prev["hdr"].get("FOCTEMP", 0)) if prev["hdr"].get("FOCPOS") == fr["hdr"].get("FOCPOS") else fr["hdr"].get("FOCTEMP", 0)
    return ("   <- " + "; ".join(w)) if w else ""


def scale_of(hdr):
    return 206.265 * float(hdr.get("XPIXSZ", 3.76)) * float(hdr.get("XBINNING", 1)) / float(hdr.get("FOCALLEN", 600))


def main():
    fix_console()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("--off", nargs="+", default=[], help="subs with compensator OFF")
    ap.add_argument("--on", nargs="+", default=[], help="subs with compensator ON")
    ap.add_argument("--map", action="store_true", help="print 3x3 field map of trailing per frame")
    args = ap.parse_args()
    entries = [(f, g) for fs, g in [(args.files, "-"), (args.off, "off"), (args.on, "on")] for f in expand(fs)]
    if not entries:
        ap.error("no files")

    frames = []
    for f, g in entries:
        print(f"measuring {os.path.basename(f)} ...", file=sys.stderr, flush=True)
        s, hdr, shape = measure(f)
        frames.append(dict(f=f, g=g, s=s, hdr=hdr, shape=shape, t=t_of(hdr), tl=t_of(hdr, "DATE-LOC"), st=shape_stats(s)))
    frames.sort(key=lambda fr: fr["t"])
    sc = scale_of(frames[0]["hdr"])

    print(f"\npixel scale {sc:.3f}\"/px.  trail = field-coherent smear inside the sub; PA in image coords (0 = +x, 90 = +y)\n")
    print(f"{'#':>2} {'grp':>3} {'local':>8} {'exp':>5} {'flt':>3} {'pier':>4} {'alt':>5} {'az':>5} {'focus':>5} {'foc C':>5} {'stars':>5} "
          + " ".join(t.rjust(n) for t, n in [('FWHM"', 6), ('minor"', 6), ('ecc', 5), ('trail px', 8), ('trail"', 6), ('"/min', 5), ('PA', 4)]))
    for i, fr in enumerate(frames):
        h, st = fr["hdr"], fr["st"]
        exp = float(h.get('EXPTIME', 0))
        print(f"{i:>2} {fr['g']:>3} {fr['tl'].strftime('%H:%M:%S'):>8} {exp:>5.0f} {h.get('FILTER', '?'):>3} "
              f"{h.get('PIERSIDE', '?')[:4]:>4} {float(h.get('CENTALT', 'nan')):>5.1f} {float(h.get('CENTAZ', 'nan')):>5.1f} "
              f"{h.get('FOCPOS', '-'):>5} {float(h.get('FOCTEMP', 'nan')):>5.1f} {st['n']:>5d} "
              f"{st['fwhm'] * sc:>6.2f} {st['fwhm_min'] * sc:>6.2f} {st['ecc']:>5.2f} {st['trail']:>8.2f} {st['trail'] * sc:>6.2f} "
              f"{st['trail'] * sc / (exp / 60):>5.2f} {st['pa']:>4.0f}" + (warn(fr, frames[i - 1] if i else None)))
        if args.map:
            for row in field_map(fr["s"], fr["shape"]):
                print("      " + "   ".join(f"{c['trail'] * sc:4.1f}\"@{c['pa']:3.0f} e{c['ecc']:.2f}" if c else "        --        " for c in row))

    if len(frames) > 1:
        print(f"\nDrift between consecutive subs (centre of field). 'pred trail' = drift rate x exposure, "
              f"compare with measured trail of the later sub.")
        print(f"{'pair':>6} {'grp':>7} {'dt min':>6} {'dx px':>7} {'dy px':>7} " + 'shift"'.rjust(7) + " " + 'rate "/min'.rjust(10) + f" {'dir':>4} "
              + " ".join(t.rjust(n) for t, n in [("rot '", 6), ('pred trail"', 11), ('meas"', 6), ('match', 5), ('rms px', 6)]))
        rows = []
        for i in range(1, len(frames)):
            a, b = frames[i - 1], frames[i]
            r = register(a["s"], b["s"], a["shape"])
            if r is None:
                print(f"{i - 1:>2}->{i:<2}  registration failed"); continue
            dt = (b["t"] - a["t"]).total_seconds() / 60
            shift = np.hypot(r["dx"], r["dy"]) * sc
            rate = shift / dt
            pred = rate * float(b["hdr"].get("EXPTIME", 0)) / 60
            pier_flip = a["hdr"].get("PIERSIDE") != b["hdr"].get("PIERSIDE")
            print(f"{i - 1:>2}->{i:<2} {a['g'] + '>' + b['g']:>7} {dt:>6.1f} {r['dx']:>7.2f} {r['dy']:>7.2f} {shift:>7.2f} {rate:>10.3f} "
                  f"{np.degrees(np.arctan2(r['dy'], r['dx'])) % 180:>4.0f} {r['rot']:>6.2f} {pred:>11.2f} {b['st']['trail'] * sc:>6.2f} "
                  f"{r['nmatch']:>5d} {r['rms']:>6.2f}" + ("   PIER FLIP" if pier_flip else ""))
            rows.append((a["g"], b["g"], rate, pier_flip, shift))
        big = [r for r in rows if r[4] > 3]
        if big:
            print(f"\n  WARNING: {len(big)} pair(s) jumped > 3\". With random directions that's dithering or a re-centre,"
                  "\n  not flexure - turn dithering off for the test or the drift numbers mean nothing.")

    # Drift fit: register every sub in a group to its first sub, fit position vs time.
    # Works with short subs - precision comes from the time span, not exposure length.
    print("\nDrift fit per group (each sub registered to the group's first sub, straight line vs time):")
    fits = {}
    for g in sorted({fr["g"] for fr in frames}):
        fs = [fr for fr in frames if fr["g"] == g]
        if len(fs) < 3:
            print(f"  {g:>3}: need >= 3 subs for a fit (have {len(fs)})"); continue
        t, xy = [], []
        for fr in fs:
            r = register(fs[0]["s"], fr["s"], fr["shape"])
            if r is None:
                continue
            mid = (fr["t"] - fs[0]["t"]).total_seconds() / 60 + float(fr["hdr"].get("EXPTIME", 0)) / 120
            t.append(mid); xy.append((r["dx"], r["dy"]))
        t, xy = np.array(t), np.array(xy) * sc
        if len(t) < 3:
            print(f"  {g:>3}: registration failed"); continue
        A = np.c_[t, np.ones_like(t)]
        (vx, x0), resx = np.linalg.lstsq(A, xy[:, 0], rcond=None)[:2]
        (vy, y0), resy = np.linalg.lstsq(A, xy[:, 1], rcond=None)[:2]
        jit = np.sqrt(np.mean(np.r_[xy[:, 0] - A @ [vx, x0], xy[:, 1] - A @ [vy, y0]] ** 2))
        err = jit / (np.sqrt(len(t)) * t.std()) if t.std() > 0 else np.nan
        rate = np.hypot(vx, vy)
        fits[g] = (rate, err)
        print(f"  {g:>3}: {len(t)} subs over {t.max() - t.min():.1f} min   drift {rate:.3f} +/- {err:.3f}\"/min   "
              f"dir {np.degrees(np.arctan2(vy, vx)) % 360:.0f} deg   jitter about the line {jit:.2f}\""
              + ("   <- jitter > 2\": dithering on, or not a steady drift" if jit > 2 else ""))
        # In-sub trail is only measurable once it clears the ~0.2 px^2 elongation floor
        # (optics + seeing) at ~3 sigma: trail >= ~3.5".
        if rate > 0:
            tmin = 3.5 / rate * 60
            print(f"       to confirm with in-sub trails at this rate: one sub of >= {tmin:.0f} s"
                  + ("  (impractical - trail this small is below the seeing floor and harmless)" if tmin > 1800 else ""))
    if "off" in fits and "on" in fits:
        (ro, eo), (rn, en) = fits["off"], fits["on"]
        print(f"\n  compensator: {ro:.3f} -> {rn:.3f}\"/min  (difference {ro - rn:+.3f} +/- {np.hypot(eo, en):.3f})")

    groups = sorted({fr["g"] for fr in frames} - {"-"})
    if groups:
        print("\nSummary by group (median):")
        for g in groups:
            fs = [fr for fr in frames if fr["g"] == g]
            rates = [r[2] for r in rows if r[0] == g and r[1] == g and not r[3]] if len(frames) > 1 else []
            print(f"  {g:>3}: {len(fs)} subs   trail {np.median([fr['st']['trail'] for fr in fs]) * sc:.2f}\" "
                  f"({np.median([fr['st']['trail'] * 60 / float(fr['hdr'].get('EXPTIME', 1)) for fr in fs]) * sc:.3f}\"/min)   "
                  f"ecc {np.median([fr['st']['ecc'] for fr in fs]):.2f}   FWHM {np.median([fr['st']['fwhm'] for fr in fs]) * sc:.2f}\"   "
                  + (f"drift {np.median(rates):.3f}\"/min ({len(rates)} pairs)" if rates else "drift n/a (need 2 consecutive subs in group)"))


if __name__ == "__main__":
    main()
