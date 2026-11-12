#!/usr/bin/env python3
"""Measure how much depth each light sub adds: sky brightness, transparency and effective time.

For every light sub (XISF or uncompressed FITS, mono or one-shot colour):

  zp        photometric zero point: the Gaia BP magnitude that gives 1 e-/s.
            The sub is plate-solved with ASTAP and its stars are matched to
            ASTAP's local Gaia star database (H18 or H17).
  sky e-    sky above the bias level per pixel, in electrons. Measured two ways:
            from the bias level (sky minus a matching master dark, else
            --pedestal, else a bias estimated from the OFFSET header for known
            cameras), and from the pixel noise (noise^2 = sky + read noise^2),
            which needs no bias level at all. An estimated bias is used only
            when the two agree within 10%; otherwise the noise value is used.
            Disagreements are reported per night.
  b/n       sky from the bias / sky from the noise. Close to 1 when gain, read
            noise and bias level are right.
  gain chk  the gain implied by sky level vs noise. It should be close to
            EGAIN; if not, the gain or the bias level is wrong.
  mu        sky brightness, mag/arcsec^2 (in that filter, on the Gaia BP scale).
  T         transparency: flux relative to the clearest subs of that filter
            (haze, thin cloud, low altitude).
  t_eff     effective time: what the sub is worth for a faint target, in
            seconds of an ideal sub (clearest sky, darkest sky of that filter,
            no read noise):  t_eff = t * T^2 * b_ref / (b + RN^2 / t)

Prints a summary per filter and night, and per filter. --csv writes every sub.

Examples:
  subdepth.py "D:/Astro/Heart Nebula" --darks "D:/Astro/Calibration Frames/Darks"
  subdepth.py "D:/Astro/Dark Shark/_selected" --gain 100=0.25 --gain 0=0.79 --csv depth.csv
  subdepth.py "D:/Astro/Old Target" --pedestal 2750 --gain 0.25
  subdepth.py lights/ --astap "C:/Program Files/astap" -v

Requirements: Python 3.9+, numpy, scipy (+ lz4/zstandard for compressed XISF),
and ASTAP (astap_cli plus the H18 or H17 star database).
"""
import argparse, csv, hashlib, json, os, re, shutil, struct, subprocess, sys, tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "subselect"))
from subselect import (read_header, read_image, find_stars, shape_stats, load_darks, match_dark,  # noqa: E402
                       temp_of, night_of, time_of, pixel_scale, fix_argv, expand, find_files)

np.seterr(all="ignore")

VERSION = "2"
CACHE = os.environ.get("SUBDEPTH_CACHE") or (
    os.path.join(os.environ.get("LOCALAPPDATA", ""), "subdepth") if os.name == "nt"
    else os.path.expanduser("~/.cache/subdepth"))
MAG_RANGE = (10.0, 15.0)    # catalogue stars used for the zero point: bright enough, rarely saturated
TILE = 128                  # sky tiles, px
SKY_PCT = 20                # sky = darkest 20% of central tiles (avoids nebula and vignetting)
CROP = 0.6                  # central fraction of the frame used for the sky
# Bias level in native ADU per OFFSET unit, measured from master darks (2026-10-31). Used when
# no master dark matches. Add a camera here once its darks have been measured.
BIAS_PER_OFFSET = (("ASI2600MM", 10.0), ("ASI533MC", 9.8))
AGREE = 1.10                # sky from bias vs from noise: warn (and distrust an estimate) beyond 10%

# ---------------------------------------------------------------- ASTAP star database (.1476)
# 36 rings of constant declination, each split into equal RA cells. File name
# <db>_<ring><cell>.1476. Each file: 110-byte text header (byte 109 = record
# size, 5 for H17/H18), then 5-byte records sorted by magnitude. A record with
# RA = 0xFFFFFF starts a magnitude group: byte 4 - 16 = magnitude * 10, byte 3
# - 128 = high byte of the declination for the records that follow.
NCELL = [1, 3, 9, 15, 21, 27, 33, 38, 43, 48, 52, 56, 60, 63, 65, 67, 68, 69,
         69, 68, 67, 65, 63, 60, 56, 52, 48, 43, 38, 33, 27, 21, 15, 9, 3, 1]
RING = 90 / 17.5


def ring_dec(r):
    """Declination range of ring r (1 = south pole, 36 = north pole), degrees."""
    if r == 1:
        return -90.0, -90 + RING / 2
    if r == 36:
        return 90 - RING / 2, 90.0
    lo = -90 + RING / 2 + (r - 2) * RING
    return lo, lo + RING


def read_tile(path):
    """Returns (ra deg, dec deg, BP mag) arrays for one .1476 file."""
    with open(path, "rb") as f:
        b = f.read()
    if len(b) < 110 or b[109] != 5:
        raise ValueError(f"{path}: not a 5-byte-record star database (H17/H18)")
    r = np.frombuffer(b, np.uint8, offset=110)[: (len(b) - 110) // 5 * 5].reshape(-1, 5).astype(np.int64)
    ra_raw = r[:, 0] | (r[:, 1] << 8) | (r[:, 2] << 16)
    hdr = ra_raw == 0xFFFFFF
    last = np.maximum.accumulate(np.where(hdr, np.arange(len(r)), -1))
    star = ~hdr & (last >= 0)
    li = last[star]
    mag = (r[li, 4] - 16) / 10.0
    dec9 = r[li, 3] - 128
    ra = ra_raw[star] * 360.0 / (2 ** 24 - 1)
    dec = ((dec9 << 16) + (r[star, 4] << 8) + r[star, 3]) * 90.0 / (2 ** 23 - 1)
    return ra, dec, mag


def find_database(astap_dir):
    for db in ("h18", "h17"):
        if os.path.exists(os.path.join(astap_dir, f"{db}_0101.1476")):
            return db
    return None


def catalog(dbdir, db, ra0, dec0, radius):
    """All catalogue stars within roughly `radius` degrees of (ra0, dec0)."""
    parts = []
    for r in range(1, 37):
        lo, hi = ring_dec(r)
        if hi < dec0 - radius or lo > dec0 + radius:
            continue
        n = NCELL[r - 1]
        top = abs(dec0) + radius
        dra = 180.0 if top > 89 else radius / np.cos(np.radians(top))
        for k in range(1, n + 1):
            a, b = (k - 1) * 360 / n, k * 360 / n
            if n == 1 or abs((ra0 - (a + b) / 2 + 180) % 360 - 180) <= (b - a) / 2 + dra:
                p = os.path.join(dbdir, f"{db}_{r:02d}{k:02d}.1476")
                if os.path.exists(p):
                    parts.append(read_tile(p))
    if not parts:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    return tuple(np.concatenate(x) for x in zip(*parts))


# ---------------------------------------------------------------- WCS (TAN, as written by ASTAP)
def world_to_pixel(w, ra, dec):
    """0-based (x = column, y = row of the array as written to FITS), plus in-front mask."""
    a0, d0 = np.radians(w["CRVAL1"]), np.radians(w["CRVAL2"])
    a, d = np.radians(ra), np.radians(dec)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.degrees(np.cos(d) * np.sin(a - a0) / cosc)
    eta = np.degrees((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc)
    inv = np.linalg.inv(np.array([[w["CD1_1"], w["CD1_2"]], [w["CD2_1"], w["CD2_2"]]]))
    x = inv[0, 0] * xi + inv[0, 1] * eta + w["CRPIX1"] - 1
    y = inv[1, 0] * xi + inv[1, 1] * eta + w["CRPIX2"] - 1
    return x, y, cosc > 0


def wcs_scale(w):
    """arcsec per pixel from the CD matrix."""
    return 3600 * np.sqrt(abs(w["CD1_1"] * w["CD2_2"] - w["CD1_2"] * w["CD2_1"]))


# ---------------------------------------------------------------- plate solving
def astap_exe(astap_dir):
    for name in ("astap_cli.exe", "astap_cli", "astap.exe", "astap"):
        p = os.path.join(astap_dir, name)
        if os.path.isfile(p):
            return p
    return shutil.which("astap_cli") or shutil.which("astap")


def default_astap_dir():
    env = os.environ.get("ASTAP_DIR")
    cands = [env] if env else []
    cands += [r"C:\Program Files\astap", "/mnt/c/Program Files/astap", "/opt/astap", "/usr/share/astap",
              "/Applications/ASTAP.app/Contents/MacOS"]
    for c in cands:
        if c and os.path.isdir(c) and find_database(c):
            return c
    return None


def _is_wsl_exe(exe):
    return os.name != "nt" and exe.lower().endswith(".exe")


def solve_dir(exe):
    """A temp folder ASTAP can read. A Windows .exe under WSL can't open \\\\wsl.localhost paths."""
    if _is_wsl_exe(exe):
        try:
            out = subprocess.run(["/mnt/c/Windows/System32/cmd.exe", "/c", "echo %TEMP%"],
                                 capture_output=True, text=True, timeout=20, cwd="/mnt/c").stdout.strip()
            d = subprocess.run(["wslpath", "-u", out], capture_output=True, text=True).stdout.strip()
            if d and os.path.isdir(d):
                return os.path.join(d, "subdepth")
        except (OSError, subprocess.SubprocessError):
            pass
    return os.path.join(tempfile.gettempdir(), "subdepth")


def write_fits(path, img):
    """Minimal float32 FITS. Array row 0 becomes FITS row 1."""
    h, w = img.shape
    cards = ["SIMPLE  =                    T", "BITPIX  =                  -32", "NAXIS   =                    2",
             f"NAXIS1  = {w:>20}", f"NAXIS2  = {h:>20}", "END"]
    head = "".join(c.ljust(80) for c in cards)
    head += " " * (-len(head) % 2880)
    data = np.ascontiguousarray(img, ">f4").tobytes()
    with open(path, "wb") as f:
        f.write(head.encode("ascii"))
        f.write(data)
        f.write(b"\0" * (-len(data) % 2880))


def _hms(v, hours):
    p = [abs(float(x)) for x in re.split(r"[ :hdms°'\"]+", v.strip()) if x not in ("", "+", "-")]
    val = p[0] + (p[1] if len(p) > 1 else 0) / 60 + (p[2] if len(p) > 2 else 0) / 3600
    val = -val if v.strip().startswith("-") else val
    return val * 15 if hours else val


def pointing(hdr):
    """Approximate (ra, dec) in degrees from the header, or None."""
    try:
        return float(hdr["RA"]), float(hdr["DEC"])
    except (KeyError, ValueError):
        pass
    try:
        return _hms(hdr["OBJCTRA"], True), _hms(hdr["OBJCTDEC"], False)
    except (KeyError, ValueError, IndexError):
        return None


def solve(img, hdr, scale, exe, tmpdir, tag):
    """Plate-solve img with ASTAP. Returns a WCS dict or None."""
    os.makedirs(tmpdir, exist_ok=True)
    base = os.path.join(tmpdir, f"s{os.getpid()}_{tag}")
    fits = base + ".fits"
    write_fits(fits, img)
    arg = fits
    if _is_wsl_exe(exe):
        arg = subprocess.run(["wslpath", "-w", fits], capture_output=True, text=True).stdout.strip()
    cmd = [exe, "-f", arg, "-z", "0"]
    if scale:
        cmd += ["-fov", f"{img.shape[0] * scale / 3600:.4f}"]
    p = pointing(hdr)
    if p:
        cmd += ["-ra", f"{p[0] / 15:.6f}", "-spd", f"{p[1] + 90:.6f}", "-r", "10"]
    else:
        cmd += ["-r", "180"]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        with open(base + ".ini", encoding="utf-8", errors="replace") as f:
            kv = dict(re.findall(r"^(\w+)=(.*?)\s*$", f.read(), re.M))
        if kv.get("PLTSOLVD", "").strip() != "T":
            return None
        return {k: float(kv[k]) for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2")}
    except (OSError, subprocess.SubprocessError, KeyError, ValueError):
        return None
    finally:
        for ext in (".fits", ".ini", ".wcs", ".log"):
            try:
                os.remove(base + ext)
            except OSError:
                pass


# ---------------------------------------------------------------- measurements
def adu_step(img):
    """Storage step of integer data: 4 when a 14-bit camera's values are stored x4 in 16 bits.
    EGAIN from the camera driver is per native ADU, so it must be divided by this."""
    v = img.ravel()[::max(1, img.size // 1_000_000)]
    if not np.all(v == np.round(v)):
        return 1
    v = v.astype(np.int64)
    v = v[v != 0]
    for k in (16, 8, 4, 2):
        if len(v) and np.all(v % k == 0):
            return k
    return 1


def superpixel(img):
    """2x2 sum: one mono pixel per Bayer cell (R+G+G+B)."""
    h, w = img.shape[0] // 2 * 2, img.shape[1] // 2 * 2
    a = img[:h, :w]
    return a[0::2, 0::2] + a[0::2, 1::2] + a[1::2, 0::2] + a[1::2, 1::2]


def sky_stats(img):
    """(sky level ADU, per-pixel noise ADU) over the darkest central tiles, stars masked."""
    h, w = img.shape
    y0, x0 = int(h * (1 - CROP) / 2), int(w * (1 - CROP) / 2)
    c = img[y0:h - y0, x0:w - x0].astype(np.float32)
    ny, nx = c.shape[0] // TILE, c.shape[1] // TILE
    if ny < 1 or nx < 1:
        raise ValueError("image too small")
    c = c[:ny * TILE, :nx * TILE]
    t = c.reshape(ny, TILE, nx, TILE)
    med = np.median(t, axis=(1, 3))
    resid = c - np.repeat(np.repeat(med, TILE, 0), TILE, 1)
    s0 = 1.4826 * np.median(np.abs(resid[::3, ::3]))
    mask = ndimage.binary_dilation(ndimage.uniform_filter(resid, 3) > 2 * s0, iterations=4)
    dx = np.full_like(c, np.nan)
    dx[:, 1:] = c[:, 1:] - c[:, :-1]
    dx[mask | np.roll(mask, 1, axis=1)] = np.nan
    lv = c.copy()
    lv[mask] = np.nan
    lvt = np.nanmedian(lv.reshape(ny, TILE, nx, TILE), axis=(1, 3))
    dxt = dx.reshape(ny, TILE, nx, TILE).transpose(0, 2, 1, 3).reshape(ny, nx, -1)
    sig = 1.4826 * np.nanmedian(np.abs(dxt - np.nanmedian(dxt, axis=2, keepdims=True)), axis=2) / np.sqrt(2)
    dark = lvt <= np.nanpercentile(lvt, SKY_PCT)
    return float(np.nanmedian(lvt[dark])), float(np.nanmedian(sig[dark]))


def photometry(img, xy, r_ap):
    """Aperture sums (ADU, local annulus background removed). NaN where it doesn't fit."""
    R = int(np.ceil(r_ap + 12)) + 1
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    h, w = img.shape
    out = np.full(len(xy), np.nan)
    for i, (sx, sy) in enumerate(xy):
        cx, cy = int(round(sx)), int(round(sy))
        if cx - R < 0 or cy - R < 0 or cx + R >= w or cy + R >= h:
            continue
        cut = img[cy - R:cy + R + 1, cx - R:cx + R + 1].astype(np.float64)
        rr = np.hypot(xx - (sx - cx), yy - (sy - cy))
        bg = np.median(cut[(rr > r_ap + 4) & (rr < r_ap + 12)])
        out[i] = (cut[rr <= r_ap] - bg).sum()
    return out


def zero_point(img, stars, w, dbdir, db, exp):
    """(zp in ADU/s, scatter, matched count) from detected stars vs catalogue."""
    h, wd = img.shape
    scale = wcs_scale(w)
    ra0, dec0 = w["CRVAL1"], w["CRVAL2"]
    radius = np.hypot(h, wd) * scale / 3600 / 2 + 0.1
    ra, dec, mag = catalog(dbdir, db, ra0, dec0, radius)
    sel = (mag >= MAG_RANGE[0]) & (mag <= MAG_RANGE[1])
    if sel.sum() < 10 or len(stars) < 10:
        return None, None, 0
    x, y, front = world_to_pixel(w, ra[sel], dec[sel])
    m = front & (x > 0) & (x < wd - 1) & (y > 0) & (y < h - 1)
    x, y, mag = x[m], y[m], mag[sel][m]
    if len(x) < 10:
        return None, None, 0
    fw = shape_stats(stars)["fwhm"]
    fw = fw if np.isfinite(fw) else 3.0
    dist, idx = cKDTree(np.c_[x, y]).query(stars[:, :2], distance_upper_bound=max(2.0, fw))
    ok = np.isfinite(dist)
    # a catalogue star claimed by two detections (close pair) is ambiguous: drop both
    uniq, cnt = np.unique(idx[ok], return_counts=True)
    ok &= ~np.isin(idx, uniq[cnt > 1])
    if ok.sum() < 10:
        return None, None, int(ok.sum())
    flux = photometry(img, stars[ok, :2], max(3.0, 1.5 * fw))
    good = flux > 0
    zp = mag[idx[ok]][good] + 2.5 * np.log10(flux[good] / exp)
    if len(zp) < 10:
        return None, None, len(zp)
    med = np.median(zp)
    mad = 1.4826 * np.median(np.abs(zp - med))
    keep = np.abs(zp - med) < 3 * max(mad, 0.02)
    return float(np.median(zp[keep])), float(1.4826 * np.median(np.abs(zp[keep] - np.median(zp[keep])))), int(keep.sum())


def measure(path, astap_dir, tmpdir):
    """Per-sub measurements in ADU (gain and darks are applied later), cached by path+size+mtime."""
    st = os.stat(path)
    key = hashlib.sha1(f"{os.path.normcase(os.path.abspath(path))}|{st.st_size}|{st.st_mtime}|{VERSION}".encode()).hexdigest()[:20]
    cpath = os.path.join(CACHE, key + ".json")
    try:
        with open(cpath, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        pass
    img, hdr = read_image(path)
    step = adu_step(img)
    bayer = bool(hdr.get("BAYERPAT", "").strip())
    if bayer:
        img = superpixel(img)
    scale = pixel_scale(hdr, None)
    if scale and bayer:
        scale *= 2
    sky_adu, sig_adu = sky_stats(img)
    stars, _ = find_stars(img / 4 if bayer else img)      # saturation limit is per raw pixel
    exe = astap_exe(astap_dir)
    db = find_database(astap_dir)
    w = solve(img, hdr, scale, exe, tmpdir, key) if exe and db else None
    zp = zp_mad = None
    n_match = 0
    if w:
        scale = wcs_scale(w)
        zp, zp_mad, n_match = zero_point(img, stars, w, astap_dir, db, float(hdr.get("EXPTIME", 0) or 0) or 1.0)
    fw = shape_stats(stars)["fwhm"] if len(stars) else np.nan
    m = dict(bayer=bayer, adu_step=step, scale=scale, solved=bool(w), zp_adu=zp, zp_mad=zp_mad, n_match=n_match,
             n_stars=len(stars), fwhm_px=float(fw) if np.isfinite(fw) else None,
             sky_adu=sky_adu, sig_adu=sig_adu, hdr=hdr)
    os.makedirs(CACHE, exist_ok=True)
    tmp = f"{cpath}.{os.getpid()}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(m, f)
    os.replace(tmp, cpath)
    return m


def _measure_safe(args):
    path = args[0]
    try:
        return path, measure(*args), None
    except Exception as e:                      # one bad file shouldn't stop the run
        return path, None, f"{type(e).__name__}: {e}"


# ---------------------------------------------------------------- electrons and depth
def estimated_pedestal(hdr, step):
    """Bias level in file ADU per raw pixel from the OFFSET header, for cameras in BIAS_PER_OFFSET; else None."""
    inst = (hdr.get("INSTRUME", "") or "").upper().replace(" ", "")
    try:
        off = float(hdr.get("OFFSET", ""))
    except ValueError:
        return None
    for model, k in BIAS_PER_OFFSET:
        if model in inst:
            return k * off * step
    return None


def agrees(ratio):
    return ratio is not None and 1 / AGREE <= ratio <= AGREE


def derive(m, exp, gain, rn, pedestal, source="dark"):
    """Electron-based values for one sub.

    pedestal: bias/dark level in ADU per output pixel, or None. source says where it came from
    ("dark", "pedestal" or "estimated bias"). The sky is measured from the pedestal when there is
    one, and always from the pixel noise too (sky_noise_e); bias_noise is their ratio.
    """
    out = dict(sky_e=None, sky_method=None, sky_noise_e=None, sky_bias_e=None, bias_noise=None,
               gain_chk=None, zp=None, rate=None, mu=None)
    if not gain or not exp:
        return out
    if m["bayer"]:
        rn = rn * 2                              # four pixels summed
    noise_e = (m["sig_adu"] * gain) ** 2 - rn ** 2
    out["sky_noise_e"] = noise_e
    out["sky_e"], out["sky_method"] = noise_e, "noise"
    if pedestal is not None:
        s = m["sky_adu"] - pedestal
        out["sky_bias_e"] = s * gain
        if noise_e > 0:
            out["bias_noise"] = s * gain / noise_e
        v = m["sig_adu"] ** 2
        if v > 0 and s > 0:
            out["gain_chk"] = (s + np.sqrt(s * s + 4 * v * rn * rn)) / (2 * v)
        out["sky_e"], out["sky_method"] = s * gain, source
    if m["zp_adu"] is not None:
        out["zp"] = m["zp_adu"] + 2.5 * np.log10(gain)
    if out["sky_e"] is not None and out["sky_e"] > 0 and m["scale"]:
        out["rate"] = out["sky_e"] / exp / m["scale"] ** 2           # e-/s/arcsec^2
        if out["zp"] is not None:
            out["mu"] = out["zp"] - 2.5 * np.log10(out["rate"])
    return out


KEEP = ("pedestal", "bias_source", "sky_bias_e", "bias_noise", "gain_chk")


def check_estimates(subs):
    """Use an estimated bias only for nights where it agrees with the noise-based sky.

    Decided per night and estimated value (median bias/noise ratio), so one night never mixes
    two methods. Subs that fall back get their noise-based values (s["_noise"]), keeping the
    bias columns for the record. Returns {(night, pedestal): median ratio} for the fallbacks.
    """
    groups = defaultdict(list)
    for s in subs:
        if s.get("bias_source") == "estimated bias" and s.get("bias_noise") is not None:
            groups[(s["night"], s["pedestal"])].append(s)
    fell = {}
    for key, g in groups.items():
        r = float(np.median([s["bias_noise"] for s in g]))
        if agrees(r):
            continue
        fell[key] = r
        for s in g:
            keep = {k: s[k] for k in KEEP if k in s}
            s.update(s["_noise"])
            s.update(keep)
    return fell


def effective_time(exp, T, rate, rate_ref, scale, rn):
    """Seconds of an ideal sub (T = 1, sky = rate_ref, no read noise) this sub is worth."""
    b = rate * scale ** 2
    return exp * T ** 2 * rate_ref * scale ** 2 / (b + rn ** 2 / exp)


def add_depth(subs, rn):
    """Fill T and t_eff, with references taken per filter from the subs themselves."""
    by_f = defaultdict(list)
    for s in subs:
        by_f[s["filter"]].append(s)
    for f, group in by_f.items():
        zps = [s["zp"] for s in group if s["zp"] is not None]
        rates = [s["rate"] for s in group if s["rate"] is not None]
        zp_ref = float(np.percentile(zps, 95)) if zps else None
        rate_ref = float(np.percentile(rates, 5)) if rates else None
        for s in group:
            s["T"] = s["t_eff"] = None
            if zp_ref is not None and s["zp"] is not None:
                s["T"] = min(1.0, 10 ** (-0.4 * (zp_ref - s["zp"])))
            if s["T"] is not None and s["rate"] and rate_ref and s["scale"]:
                r = rn * 2 if s["bayer"] else rn
                s["t_eff"] = min(s["exp"], effective_time(s["exp"], s["T"], max(s["rate"], rate_ref),
                                                          rate_ref, s["scale"], r))


def fmt(v, w, d=1):
    return f"{v:{w}.{d}f}" if isinstance(v, (int, float)) and v is not None and np.isfinite(v) else " " * (w - 1) + "-"


def summarise(rows, keyf, label):
    groups = defaultdict(list)
    for s in rows:
        groups[keyf(s)].append(s)
    print(f"\n{label:<24} {'subs':>5} {'hours':>6} {'eff h':>6} {'eff %':>5} {'mu':>6} {'T':>5} {'sky e-':>7} {'zp':>6}")
    for k in sorted(groups, key=lambda k: tuple(str(x) for x in k)):
        g = groups[k]
        hours = sum(s["exp"] for s in g) / 3600
        eff = [s["t_eff"] for s in g if s["t_eff"] is not None]
        eh = sum(eff) / 3600 if eff else None
        med = lambda key: float(np.median([s[key] for s in g if s[key] is not None])) if any(s[key] is not None for s in g) else None
        name = " ".join(str(x) for x in k)
        pct = 100 * eh / hours if eh is not None and hours else None
        print(f"{name:<24} {len(g):>5} {hours:6.2f} {fmt(eh, 6, 2)} {fmt(pct, 5, 0)} {fmt(med('mu'), 6, 2)} "
              f"{fmt(med('T'), 5, 2)} {fmt(med('sky_e'), 7, 1)} {fmt(med('zp'), 6, 2)}")


def gain_for(specs, hdr):
    """e-/ADU from --gain values: a bare number, or SETTING=VALUE matched to the GAIN header."""
    plain = None
    for spec in specs:
        k, eq, v = spec.partition("=")
        if not eq:
            plain = float(k)
            continue
        try:
            if float(k) == float(hdr.get("GAIN", "nan")):
                return float(v)
        except ValueError:
            pass
    return plain


def main(argv=None):
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="light files or folders (searched recursively)")
    ap.add_argument("--darks", action="append", default=[],
                    help="folder of master darks (matched by exposure and temperature); can be given more than once. "
                         "Without darks, the sky comes from the pixel noise")
    ap.add_argument("--pedestal", type=float, metavar="ADU",
                    help="bias level in ADU (per pixel, as in the file) for subs without a matching dark. "
                         "Default: estimated from the OFFSET header for known cameras")
    ap.add_argument("--read-noise", type=float, default=1.4, help="camera read noise, e- (default %(default)s)")
    ap.add_argument("--gain", action="append", default=[], metavar="[SETTING=]E_PER_ADU",
                    help="e-/ADU for frames without EGAIN. Either one value, or one per camera gain setting "
                         "(the GAIN header), e.g. --gain 100=0.25 --gain 0=0.79")
    ap.add_argument("--astap", help="ASTAP folder (astap_cli and the H18/H17 star database). "
                                    "Default: $ASTAP_DIR, then C:\\Program Files\\astap and /opt/astap")
    ap.add_argument("--csv", help="write the per-sub table to this CSV file")
    ap.add_argument("-v", "--verbose", action="store_true", help="print one line per sub")
    ap.add_argument("-j", "--jobs", type=int, default=min(4, os.cpu_count() or 1),
                    help="parallel workers (default %(default)s)")
    a = ap.parse_args(fix_argv(sys.argv[1:] if argv is None else list(argv)))
    a.paths, a.darks = expand(a.paths), expand(a.darks)
    a.jobs = max(1, min(a.jobs, 61) if os.name == "nt" else a.jobs)
    astap_dir = a.astap or default_astap_dir()
    if not astap_dir or not find_database(astap_dir) or not astap_exe(astap_dir):
        sys.exit("ASTAP with the H18 or H17 star database not found; give its folder with --astap")
    tmpdir = solve_dir(astap_exe(astap_dir))

    files = []
    for p in find_files(a.paths, set()):
        try:
            h = read_header(p)
        except Exception as e:
            print(f"  skipping {p}: {e}", file=sys.stderr)
            continue
        if "LIGHT" in h.get("IMAGETYP", "LIGHT").upper() and "MASTER" not in h.get("IMAGETYP", "").upper():
            files.append(p)
    if not files:
        sys.exit("no light frames found")
    print(f"{len(files)} light frames; ASTAP database {find_database(astap_dir).upper()} in {astap_dir}")

    darks = load_darks(a.darks) if a.darks else []
    subs, done, scaled = [], 0, set()
    jobs = [(p, astap_dir, tmpdir) for p in files]
    ex = ProcessPoolExecutor(a.jobs) if a.jobs > 1 else None
    try:
        for path, m, err in (ex.map(_measure_safe, jobs) if ex else map(_measure_safe, jobs)):
            done += 1
            if done % 25 == 0 or done == len(files):
                print(f"  measured {done}/{len(files)}", file=sys.stderr)
            if err:
                print(f"  {path}: {err}", file=sys.stderr)
                continue
            hdr = m["hdr"]
            exp = float(hdr.get("EXPTIME", 0) or hdr.get("EXPOSURE", 0) or 0)
            egain = float(hdr.get("EGAIN", 0) or 0) / m["adu_step"]
            gain = egain or gain_for(a.gain, hdr)
            if egain and m["adu_step"] > 1:
                scaled.add(m["adu_step"])
            ped, source = None, None
            if darks and exp:
                ped = match_dark(darks, exp, temp_of(hdr))
                source = "dark" if ped is not None else None
            if ped is None and a.pedestal is not None:
                ped, source = a.pedestal, "pedestal"
            if ped is None:
                ped = estimated_pedestal(hdr, m["adu_step"])
                source = "estimated bias" if ped is not None else None
            if ped is not None and m["bayer"]:
                ped *= 4                                 # superpixel = four pixels summed
            s = dict(file=os.path.basename(path), path=path, night=night_of(path, hdr), time=time_of(hdr),
                     filter=(hdr.get("FILTER", "") or "none").strip() or "none", exp=exp, bayer=m["bayer"],
                     scale=m["scale"], solved=m["solved"], n_match=m["n_match"], zp_mad=m["zp_mad"],
                     fwhm_as=m["fwhm_px"] * m["scale"] if m["fwhm_px"] and m["scale"] else None,
                     alt=float(hdr["CENTALT"]) if hdr.get("CENTALT") not in (None, "") else None)
            s["pedestal"] = ped
            s.update(derive(m, exp, gain, a.read_noise, ped, source or "dark"))
            s["bias_source"] = source
            if source == "estimated bias":
                s["_noise"] = derive(m, exp, gain, a.read_noise, None)
            subs.append(s)
    finally:
        if ex:
            ex.shutdown()
    if not subs:
        sys.exit("nothing measured")
    for k in sorted(scaled):
        print(f"  EGAIN divided by {k}: pixel values are stored in steps of {k} (e.g. 14-bit camera in 16-bit files)")
    nogain = [s for s in subs if s["sky_method"] is None]
    if nogain:
        print(f"  {len(nogain)} subs have no EGAIN and no matching --gain: no electrons, no t_eff", file=sys.stderr)
    check_estimates(subs)
    for s in subs:
        s["eff"] = (s["sky_e"] / (s["sky_e"] + (a.read_noise * (2 if s["bayer"] else 1)) ** 2) * 100
                    if s["sky_e"] and s["sky_e"] > 0 else None)
    add_depth(subs, a.read_noise)
    subs.sort(key=lambda s: (s["night"], s["time"] is None, s["time"] or datetime.min, s["file"]))

    unsolved = [s for s in subs if not s["solved"]]
    if unsolved:
        print(f"  {len(unsolved)} subs did not plate-solve (no zero point, no t_eff): "
              + ", ".join(s["file"] for s in unsolved[:5]) + (" ..." if len(unsolved) > 5 else ""))
    counts = defaultdict(int)
    for s in subs:
        if s["sky_method"]:
            counts[s["sky_method"]] += 1
    names = {"dark": "master dark", "pedestal": "--pedestal", "estimated bias": "bias estimated from OFFSET",
             "noise": "pixel noise"}
    print("  sky from: " + ", ".join(f"{names[k]} {n}" for k, n in sorted(counts.items(), key=lambda kv: -kv[1])))
    ratios = defaultdict(list)
    for s in subs:
        if s["bias_noise"] is not None:
            ratios[(s["night"], s["bias_source"])].append(s["bias_noise"])
    bad = {k: float(np.median(r)) for k, r in ratios.items() if not agrees(float(np.median(r)))}
    if bad:
        print(f"  sky from the bias and from the pixel noise disagree by more than {(AGREE - 1) * 100:.0f}% "
              "(check the gain, read noise, offset or darks):")
        for (night, src), r in sorted(bad.items()):
            used = "used the noise value" if src == "estimated bias" else f"used the {names[src]}"
            print(f"    {night}: {len(ratios[(night, src)])} subs, bias/noise {r:.2f} ({names[src]}; {used})")

    if a.verbose:
        print(f"\n{'night':<10} {'time':<8} {'filter':<6} {'exp':>5} {'alt':>4} {'FWHM':>5} {'stars':>5} {'zp':>6} "
              f"{'sky e-':>7} {'b/n':>4} {'gain?':>5} {'mu':>6} {'T':>5} {'eff':>4} {'t_eff':>6}  file")
        for s in subs:
            tm = s["time"].strftime("%H:%M:%S") if s["time"] else ""
            print(f"{s['night']:<10} {tm:<8} {s['filter']:<6} {s['exp']:5.0f} {fmt(s['alt'], 4, 0)} {fmt(s['fwhm_as'], 5, 2)} "
                  f"{s['n_match']:5d} {fmt(s['zp'], 6, 2)} {fmt(s['sky_e'], 7, 1)} {fmt(s['bias_noise'], 4, 2)} "
                  f"{fmt(s['gain_chk'], 5, 3)} "
                  f"{fmt(s['mu'], 6, 2)} {fmt(s['T'], 5, 2)} {fmt(s['eff'], 3, 0)}% {fmt(s['t_eff'], 6, 0)}  {s['file']}")

    summarise(subs, lambda s: (s["filter"], s["night"]), "filter night")
    summarise(subs, lambda s: (s["filter"],), "filter (all nights)")
    print("\nmu = sky mag/arcsec^2 (Gaia BP scale), T = transparency vs clearest subs, "
          "eff h = hours of ideal subs (clearest, darkest sky of that filter, no read noise)")

    if a.csv:
        cols = ["file", "night", "time", "filter", "exp", "alt", "scale", "fwhm_as", "solved", "n_match", "zp",
                "zp_mad", "sky_e", "sky_method", "pedestal", "sky_bias_e", "sky_noise_e", "bias_noise", "gain_chk",
                "rate", "mu", "T", "eff", "t_eff", "path"]
        with open(a.csv, "w", newline="", encoding="utf-8-sig") as f:
            wr = csv.writer(f)
            wr.writerow(cols)
            for s in subs:
                wr.writerow([(f"{s[c]:.5g}" if isinstance(s[c], float) else s[c]) if s[c] is not None else "" for c in cols])
        print(f"wrote {a.csv}")


if __name__ == "__main__":
    main()
