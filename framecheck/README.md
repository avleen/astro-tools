<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# framecheck

Checks that light subs are framed where they should be. Each sub is
plate-solved with ASTAP, and its centre and rotation are compared with the
target the sequence asked for.

Use it after a night with centring problems, such as a failed plate solve or a
sync that went wrong, to see whether the subs are still usable. It also checks
mosaic panels before stacking.

## Install

Same as [subdepth](../subdepth): numpy, scipy, lz4, zstandard, plus
[ASTAP](https://www.hnsky.org/astap.htm) with the H18 or H17 star database.
ASTAP is looked for in the same places.

## Use

```
python framecheck.py "D:\Astro\Veil Panel 1\NIGHT_2026-11-11"
python framecheck.py "D:\Astro\Veil Panel 1\NIGHT_2026-11-11" -v
python framecheck.py "D:\Astro\Heart Nebula" --max-shift 3
python framecheck.py lights\ --ra "20 49 32" --dec "+29 23 25" --pa 165
```

- **Where the target comes from:** each sub's header, as N.I.N.A. writes it.
  That's `OBJCTRA` and `OBJCTDEC` for the target centre and `OBJCTROT` for the
  position angle. `--ra`, `--dec` and `--pa` override them for every file.
- **Speed:** subs are binned 2×2 before solving. 28 ASI2600MM subs took
  14 s with the default 4 workers.
- **What it prints:** one summary line per target (`OBJECT`), plus one line
  for each flagged sub. With `-v` it prints a line for every sub.
- **Exit code:** 1 if any sub is flagged, so it can be used in a script.

```
Veil Panel 1: 10 subs, first 0.4', median 0.5', max 0.5', rotation within 0.1 deg; 0 flagged
  time                 flt  offset    x %    y %     PA    dPA
  2026-11-11 20:03:03    G    0.4'   -0.2   +0.3  164.9   -0.1
```

| column | meaning |
|---|---|
| offset | distance of the frame centre from the target, arcmin |
| x %, y % | the same offset as a share of the frame width and height |
| PA | solved position angle, in N.I.N.A.'s `OBJCTROT` convention |
| dPA | difference from the target PA, modulo 180°. A frame turned 180° covers the same sky. |

A sub is flagged when x or y is beyond `--max-shift` (default 5%), or dPA is
beyond `--max-rot` (default 1°).

Dithering moves each sub by a fraction of a percent, and the moves add up over
a night. Drift of 1–2% of the frame is normal, and stacking registration
handles it at the cost of a thin crop. A sub flagged at 5% or more loses that
much of the frame, and on a mosaic it eats into the overlap with the next
panel.
