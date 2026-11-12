#!/usr/bin/env python3
"""Check that light subs are framed where they should be.

Plate-solves each sub with ASTAP and compares its centre and rotation with the
target the sequence asked for. Use it after a night with centring problems
(a failed plate solve, a sync that went wrong) to see whether the subs are
still usable, or to check mosaic panels before stacking.

The target comes from each sub's header, as N.I.N.A. writes it:
OBJCTRA / OBJCTDEC (target centre) and OBJCTROT (position angle). Give --ra,
--dec and --pa to use other values for every file.

Per sub it prints:

  offset     distance of the frame centre from the target, arcmin
  x, y       the same offset as a share of the frame width / height, %
             (signed: + is towards higher column / row numbers)
  PA, dPA    solved position angle (N.I.N.A.'s OBJCTROT convention) and the
             difference from the target. A frame turned 180 degrees covers the
             same sky, so dPA is reported modulo 180.

Subs are flagged when x or y is beyond --max-shift (default 5%) or dPA beyond
--max-rot (default 1 degree). Dithering moves subs by a fraction of a percent
each, so drift of 1-2% over a night is normal and costs only a little crop.

Examples:
  framecheck.py "D:\\Astro\\Veil Panel 1\\NIGHT_2026-11-11"
  framecheck.py "D:\\Astro\\Heart Nebula" --max-shift 3
  framecheck.py lights\\ --ra "20 49 32" --dec "+29 23 25" --pa 165 -v

Requirements: Python 3.9+, numpy, scipy, lz4/zstandard for compressed XISF,
and ASTAP with the H18 or H17 star database (found as for subdepth).
"""
import argparse, math, os, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "subselect"))
sys.path.insert(0, os.path.join(HERE, "..", "subdepth"))
import numpy as np
from subselect import read_header, read_image, find_files, expand, fix_argv, pixel_scale, time_of
from subdepth import (solve, solve_dir, astap_exe, find_database, default_astap_dir,
                      world_to_pixel, wcs_scale, _hms)


def target_of(hdr, ra=None, dec=None, pa=None):
    """(ra_deg, dec_deg, pa_deg or None) from the options, else the header; None if unknown."""
    try:
        r = _hms(ra, True) if ra else _hms(str(hdr["OBJCTRA"]).strip("' "), True)
        d = _hms(dec, False) if dec else _hms(str(hdr["OBJCTDEC"]).strip("' "), False)
    except (KeyError, ValueError, IndexError):
        return None
    p = pa if pa is not None else hdr.get("OBJCTROT")
    try:
        p = float(p) if p not in (None, "") else None
    except ValueError:
        p = None
    return r, d, p


def position_angle(w):
    """Position angle of the frame, degrees, in N.I.N.A.'s OBJCTROT convention."""
    return math.degrees(math.atan2(w["CD1_2"], w["CD2_2"])) % 360


def framing(w, shape, target):
    """Offset of the frame centre from the target: arcmin, x %, y %, PA, dPA (mod 180)."""
    h, wd = shape
    tx, ty, front = world_to_pixel(w, np.array([target[0]]), np.array([target[1]]))
    if not front[0]:
        return None
    dx, dy = (wd - 1) / 2 - tx[0], (h - 1) / 2 - ty[0]
    pa = position_angle(w)
    dpa = None
    if target[2] is not None:
        dpa = (pa - target[2] + 90) % 180 - 90
    return dict(offset=math.hypot(dx, dy) * wcs_scale(w) / 60, x=100 * dx / wd, y=100 * dy / h, pa=pa, dpa=dpa)


def bin2(img):
    h, w = img.shape
    return img[: h // 2 * 2, : w // 2 * 2].reshape(h // 2, 2, w // 2, 2).mean(axis=(1, 3))


def check(path, astap_dir, tmpdir, opts):
    img, hdr = read_image(path)
    tgt = target_of(hdr, *opts)
    if tgt is None:
        return dict(error="no target in header (OBJCTRA/OBJCTDEC); give --ra and --dec")
    scale = pixel_scale(hdr, None)
    small = bin2(img)                           # 4x fewer pixels, solves just as well
    tag = f"{abs(hash(path)) % 10**8}"
    w = solve(small, hdr, scale * 2 if scale else None, astap_exe(astap_dir), tmpdir, tag)
    if not w:
        return dict(error="not solved")
    f = framing(w, small.shape, tgt)
    if f is None:
        return dict(error="target not in front of the frame")
    t = time_of(hdr)
    return dict(f, object=str(hdr.get("OBJECT", "")).strip("' ") or "(no OBJECT)",
                filter=str(hdr.get("FILTER", "")).strip("' ") or "-",
                time=t.isoformat(sep=" ") if t else "", file=os.path.basename(path), target_pa=tgt[2])


def _check_safe(args):
    try:
        return args[0], check(*args)
    except Exception as e:                      # one bad file shouldn't stop the run
        return args[0], dict(error=f"{type(e).__name__}: {e}")


def main(argv=None):
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="light files or folders (searched recursively)")
    ap.add_argument("--ra", help='target RA for every file, hours ("20 49 32" or 20:49:32)')
    ap.add_argument("--dec", help='target Dec for every file, degrees ("+29 23 25")')
    ap.add_argument("--pa", type=float, help="target position angle for every file, degrees")
    ap.add_argument("--max-shift", type=float, default=5.0, help="flag subs off target by more than this %% of the frame (default %(default)s)")
    ap.add_argument("--max-rot", type=float, default=1.0, help="flag subs rotated by more than this, degrees (default %(default)s)")
    ap.add_argument("--astap", help="ASTAP folder. Default: $ASTAP_DIR, then C:\\Program Files\\astap and /opt/astap")
    ap.add_argument("-v", "--verbose", action="store_true", help="print one line per sub (default: flagged subs only)")
    ap.add_argument("-j", "--jobs", type=int, default=min(4, os.cpu_count() or 1), help="parallel workers (default %(default)s)")
    a = ap.parse_args(fix_argv(sys.argv[1:] if argv is None else list(argv)))
    if (a.ra is None) != (a.dec is None):
        sys.exit("give both --ra and --dec, or neither")
    astap_dir = a.astap or default_astap_dir()
    if not astap_dir or not find_database(astap_dir) or not astap_exe(astap_dir):
        sys.exit("ASTAP with the H18 or H17 star database not found; give its folder with --astap")
    tmpdir = solve_dir(astap_exe(astap_dir))

    files = []
    for p in find_files(expand(a.paths), set()):
        try:
            h = read_header(p)
        except Exception as e:
            print(f"  skipping {p}: {e}", file=sys.stderr)
            continue
        if "LIGHT" in h.get("IMAGETYP", "LIGHT").upper() and "MASTER" not in h.get("IMAGETYP", "").upper():
            files.append(p)
    if not files:
        sys.exit("no light frames found")
    print(f"{len(files)} light frames")

    opts = (a.ra, a.dec, a.pa)
    jobs = [(p, astap_dir, tmpdir, opts) for p in files]
    ex = ProcessPoolExecutor(a.jobs) if a.jobs > 1 and len(files) > 1 else None
    try:
        results = list(ex.map(_check_safe, jobs) if ex else map(_check_safe, jobs))
    finally:
        if ex:
            ex.shutdown()

    groups, failed = defaultdict(list), []
    for path, r in results:
        if "error" in r:
            failed.append((path, r["error"]))
        else:
            groups[r["object"]].append(r)

    head = f"  {'time':19} {'flt':>4} {'offset':>7} {'x %':>6} {'y %':>6} {'PA':>6} {'dPA':>6}"
    flagged_total = 0
    for obj in sorted(groups):
        subs = sorted(groups[obj], key=lambda s: s["time"])
        for s in subs:
            s["flag"] = abs(s["x"]) > a.max_shift or abs(s["y"]) > a.max_shift or \
                (s["dpa"] is not None and abs(s["dpa"]) > a.max_rot)
        flagged = [s for s in subs if s["flag"]]
        flagged_total += len(flagged)
        off = [s["offset"] for s in subs]
        rot = [abs(s["dpa"]) for s in subs if s["dpa"] is not None]
        print(f"\n{obj}: {len(subs)} subs, first {off[0]:.1f}', median {np.median(off):.1f}', max {max(off):.1f}'"
              + (f", rotation within {max(rot):.1f} deg" if rot else ", no target PA") + f"; {len(flagged)} flagged")
        show = subs if a.verbose else flagged
        if show:
            print(head)
            for s in show:
                dpa = f"{s['dpa']:+6.1f}" if s["dpa"] is not None else f"{'-':>6}"
                print(f"  {s['time'][:19]:19} {s['filter']:>4} {s['offset']:6.1f}' {s['x']:+6.1f} {s['y']:+6.1f} "
                      f"{s['pa']:6.1f} {dpa}" + ("  <-- off target" if s["flag"] else ""))
    for path, err in failed:
        print(f"  {path}: {err}", file=sys.stderr)
    print(f"\n{sum(len(g) for g in groups.values())} solved, {len(failed)} not, {flagged_total} flagged "
          f"(beyond {a.max_shift:g}% of the frame or {a.max_rot:g} deg)")
    return 1 if flagged_total else 0


if __name__ == "__main__":
    sys.exit(main())
