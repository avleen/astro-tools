#!/usr/bin/env python3
"""Recover images from PixInsight swap files after a crash.

PixInsight keeps each open image's undo history as swap files in its swap
directories (Edit > Global Preferences > Directories and Network; by default
the system temp folder, %TEMP% on Windows), named

  ~PI~<16-char id>~<view id>.swp            128-byte header + first chunk of pixels
  ~PI~<16-char id>~<view id>-partNNN.swp    further chunks (one per swap directory)

They are normally deleted on exit, but survive a crash. Each set is one history
state: the image as it was *before* a process was applied to it. So the newest
state of a view is one step behind what was on screen when PixInsight died.

This converts every complete state to an uncompressed XISF that PixInsight opens
normally, named <time>_<view id>_<id>.xisf, plus an optional auto-stretched PNG
preview so you can see which is which.

COPY THE SWAP FILES SOMEWHERE SAFE FIRST (or use --copy-to) - a running
PixInsight may delete them.

Usage:
  pi_swap_recover.py [SWAP_DIR ...] OUT_DIR [--previews] [--view NAME] [--latest] [--copy-to DIR]

With no SWAP_DIR, the swap directories from PixInsight's settings are used
(or the system temp folder if none are set).
"""
import argparse, glob, os, re, shutil, struct, sys, tempfile
from datetime import datetime
from xml.sax.saxutils import quoteattr

import numpy as np

HDR = 128
MAGIC = b"egamiwar"      # "rawimage" reversed
SWAP_GLOB = "~PI~*~*.swp"


def parse_header(b):
    if b[:8] != MAGIC:
        raise ValueError("not a PixInsight swap file")
    is_float, bits = b[12], b[15]
    w, h, nc, cs = struct.unpack_from("<4I", b, 0x1C)
    parts = struct.unpack_from("<I", b, 0x68)[0]
    dtype = {(1, 32): "<f4", (1, 64): "<f8", (0, 8): "u1", (0, 16): "<u2", (0, 32): "<u4"}.get((is_float, bits))
    if dtype is None:
        raise ValueError(f"unknown sample format float={is_float} bits={bits}")
    return dict(w=w, h=h, nc=nc, rgb=cs == 1, parts=max(parts, 1), dtype=dtype)


def pi_swap_dirs():
    """Swap directories from PixInsight's settings, else the system temp dir (PixInsight's default)."""
    base = os.environ.get("APPDATA") if os.name == "nt" else os.path.expanduser("~/.config")
    dirs = []
    try:
        with open(os.path.join(base or "", "Pleiades", "PixInsight.ini"), encoding="utf-8", errors="replace") as f:
            for m in re.finditer(r"^Global\\Preferences\\ImageWindow\\SwapDirectory\d+=(.+?)\s*$", f.read(), re.M):
                dirs.append(m.group(1).strip('"'))
    except OSError:
        pass
    return dirs or [tempfile.gettempdir()]


def clean_dirs(dirs):
    out = []                                     # "D:\Temp\" in cmd arrives as D:\Temp" - drop the stray quote
    for d in dirs:
        d = os.path.abspath(os.path.expanduser(os.path.expandvars(d.strip().strip('"'))))
        if os.path.normcase(d) not in map(os.path.normcase, out):
            out.append(d)
    return out


def swap_files(dirs):
    return sorted({p for d in dirs for p in glob.glob(os.path.join(glob.escape(d), SWAP_GLOB))})


def find_states(swap_dirs):
    if isinstance(swap_dirs, (str, os.PathLike)):
        swap_dirs = [swap_dirs]
    states = {}
    for p in swap_files(swap_dirs):              # parts may be spread over several swap directories
        m = re.match(r"~PI~([A-Z0-9]+)~(.+?)(?:-part(\d+))?\.swp$", os.path.basename(p), re.I)
        if not m:
            continue
        sid, view, part = m.group(1), m.group(2), int(m.group(3) or 0)
        states.setdefault(sid, dict(id=sid, view=view, files={}))["files"][part] = p
    for s in states.values():
        f = s["files"].get(0) or max(s["files"].values(), key=os.path.getmtime)
        s["mtime"] = os.path.getmtime(f)
    return sorted(states.values(), key=lambda s: s["mtime"])


def load_state(s):
    if 0 not in s["files"]:
        raise ValueError("main file (without -partNNN) missing")
    with open(s["files"][0], "rb") as f:
        hdr = parse_header(f.read(HDR))
        chunks = [f.read()]
    missing = [i for i in range(1, hdr["parts"]) if i not in s["files"]]
    if missing:
        raise ValueError(f"missing part(s) {missing}")
    for i in range(1, hdr["parts"]):
        with open(s["files"][i], "rb") as f:
            chunks.append(f.read())
    data = np.frombuffer(b"".join(chunks), hdr["dtype"])
    n = hdr["w"] * hdr["h"] * hdr["nc"]
    if data.size < n:
        raise ValueError(f"truncated: {data.size} of {n} samples (swap compression on?)")
    # The pixels are split evenly across the parts, and each part holds its share
    # of every channel in turn: part k = [ch0 slice k, ch1 slice k, ...]
    parts, nc = hdr["parts"], hdr["nc"]
    if n % (parts * nc):
        raise ValueError(f"{n} samples don't split evenly over {parts} parts")
    data = data[:n].reshape(parts, nc, n // parts // nc).transpose(1, 0, 2)
    return data.reshape(nc, hdr["h"], hdr["w"]), hdr


def write_xisf(path, img, rgb, note):
    nc, h, w = img.shape
    fmt = ("Float" if img.dtype.kind == "f" else "UInt") + str(img.dtype.itemsize * 8)
    bounds = ' bounds="0:1"' if fmt.startswith("Float") else ""
    data = np.ascontiguousarray(img, img.dtype.newbyteorder("<")).tobytes()
    off = 65536                                  # fixed, so the XML can state it before it's written
    xml = ('<?xml version="1.0" encoding="UTF-8"?>\n'
           '<xisf version="1.0" xmlns="http://www.pixinsight.com/xisf" '
           'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
           'xsi:schemaLocation="http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd">\n'
           f'<Image geometry="{w}:{h}:{nc}" sampleFormat="{fmt}"{bounds} '
           f'colorSpace="{"RGB" if rgb else "Gray"}" location="attachment:{off}:{len(data)}">\n'
           f'<FITSKeyword name="HISTORY" value="" comment={quoteattr(note)}/>\n'
           '</Image>\n'
           '<Metadata><Property id="XISF:CreatorApplication" type="String">pi_swap_recover</Property></Metadata>\n'
           '</xisf>\n').encode("utf-8")
    head = b"XISF0100" + struct.pack("<II", len(xml), 0) + xml
    if len(head) > off:
        raise ValueError("XML header too long")
    with open(path, "wb") as f:
        f.write(head + b"\0" * (off - len(head)))
        f.write(data)


def preview(path, img):
    from PIL import Image
    x = img[:, ::4, ::4].astype(np.float32)
    if img.dtype.kind in "ui":                   # integer images: scale to 0..1 first
        x /= np.iinfo(img.dtype).max
    out = []
    for c in x:                                  # per-channel STF-style auto stretch
        med = np.median(c)
        mad = 1.4826 * np.median(np.abs(c - med)) or 1e-6
        lo = min(max(med - 2.8 * mad, 0), 0.999)
        c = np.clip((c - lo) / (1 - lo), 0, 1)
        m = np.median(c)
        t = 0.25
        mtf = (t - 1) * m / ((2 * t - 1) * m - t) if 0 < m < 1 else 0.5   # puts the median at t
        out.append(((mtf - 1) * c) / ((2 * mtf - 1) * c - mtf))
    a = np.clip(np.stack(out, -1) * 255, 0, 255).astype(np.uint8)
    Image.fromarray(a[..., 0] if a.shape[-1] == 1 else a).save(path)


def safe_name(s, maxlen=64):
    """Usable as a file name on Windows too (view ids normally are, but be sure)."""
    s = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", s)[:maxlen].rstrip(" .")
    return s or "_"


def long_path(p):
    p = os.path.abspath(p)                       # get past MAX_PATH (260) on Windows
    return "\\\\?\\" + p if os.name == "nt" and len(p) >= 250 and not p.startswith("\\\\") else p


def copy_swap(dirs, dest):
    os.makedirs(dest, exist_ok=True)
    n = 0
    for p in swap_files(dirs):
        if os.path.normcase(os.path.dirname(p)) == os.path.normcase(os.path.abspath(dest)):
            continue
        try:
            shutil.copy2(p, os.path.join(dest, os.path.basename(p)))   # copy2 keeps the mtime, which orders states
            n += 1
        except PermissionError:
            print(f"NOT COPIED (locked - is PixInsight still running?): {p}")
    print(f"copied {n} swap file(s) to {dest}")


def main(argv=None):
    for st in (sys.stdout, sys.stderr):          # never die on a cp1252 console or redirected output
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("swap_dir", nargs="*", help="swap directories (default: from PixInsight's settings)")
    ap.add_argument("out_dir")
    ap.add_argument("--previews", action="store_true", help="also write auto-stretched PNG previews (needs pillow)")
    ap.add_argument("--view", help="only this view id (e.g. MasterLight_H)")
    ap.add_argument("--latest", action="store_true", help="only the newest state of each view")
    ap.add_argument("--copy-to", metavar="DIR", help="first copy the swap files here (keeping mtimes), then recover from the copies")
    a = ap.parse_intermixed_args(argv)

    dirs = clean_dirs(a.swap_dir or pi_swap_dirs())
    if not a.swap_dir:
        print("swap directories (from PixInsight's settings):", ", ".join(dirs))
    if a.copy_to:
        copy_swap(dirs, a.copy_to)
        dirs = clean_dirs([a.copy_to])
    states = find_states(dirs)
    if not states:
        print(f"no {SWAP_GLOB} files in {', '.join(dirs)}\n"
              "PixInsight's swap directories are set in Edit > Global Preferences > Directories and Network "
              "(default: the system temp folder, %TEMP% on Windows).")
        return 1
    if a.view:
        states = [s for s in states if s["view"] == a.view]
        if not states:
            print(f"no states for view {a.view!r}")
            return 1
    if a.latest:
        newest = {}
        for s in states:
            newest[s["view"]] = s
        states = sorted(newest.values(), key=lambda s: s["mtime"])
    os.makedirs(a.out_dir, exist_ok=True)
    print(f"{'saved at':19}  {'view':24} {'size':>15}  {'min':>8} {'median':>8} {'max':>8}  file")
    for s in states:
        t = datetime.fromtimestamp(s["mtime"])
        try:
            img, hdr = load_state(s)
        except PermissionError as e:
            print(f"{t:%Y-%m-%d %H:%M:%S}  {s['view']:24} SKIPPED: locked, is PixInsight still running? "
                  f"Close it or copy the swap files first ({e.filename})")
            continue
        except Exception as e:
            print(f"{t:%Y-%m-%d %H:%M:%S}  {s['view']:24} SKIPPED: {e}")
            continue
        name = f"{t:%Y%m%d_%H%M%S}_{safe_name(s['view'])}_{s['id']}"
        base = long_path(os.path.join(a.out_dir, name))
        try:
            write_xisf(base + ".xisf", img, hdr["rgb"],
                       f"recovered from PixInsight swap state {s['id']} of {s['view']} saved {t:%Y-%m-%d %H:%M:%S}")
            if a.previews:
                preview(base + ".png", img)
        except OSError as e:
            print(f"{t:%Y-%m-%d %H:%M:%S}  {s['view']:24} NOT WRITTEN: {e}")
            continue
        smp = img[:, ::8, ::8]
        print(f"{t:%Y-%m-%d %H:%M:%S}  {s['view']:24} {hdr['w']:>5}x{hdr['h']}x{hdr['nc']}  "
              f"{smp.min():8.5f} {np.median(smp):8.5f} {smp.max():8.5f}  {name}.xisf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
