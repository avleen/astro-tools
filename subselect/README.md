<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# subselect

Grades astrophotography light subs (XISF from N.I.N.A., or uncompressed FITS) on
star size, shape, star count and background. It can then copy the keepers into a
folder tree that PixInsight WBPP can pair with same-night flats.

## Install

Linux / macOS / WSL:

```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Windows (cmd.exe or PowerShell, Python from python.org):

```
py -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Use

```
python subselect.py "D:/Astro/Heart Nebula"                 # report only, changes nothing
python subselect.py "D:/Astro/Heart Nebula" --copy          # copy to "D:/Astro/Heart Nebula/_selected"
python subselect.py lights/ --darks masters/ --csv out.csv  # add sky e- / read-noise efficiency, CSV
```

On Windows:

```
python subselect.py "D:\Astro\Heart Nebula"
python subselect.py "D:\Astro\Heart Nebula" --copy --out "E:\Stacking\Heart"
python subselect.py "D:\Astro\Heart Nebula\NIGHT_*\LIGHT" --out "D:\Astro\Heart Nebula\_selected"
python subselect.py D:\Astro\Heart\*.xisf --csv report.csv
```

- Put paths that contain spaces in double quotes. Forward slashes work too.
- cmd.exe and PowerShell do not expand wildcards, so `subselect.py` expands
  `*`, `?` and `[...]` itself (and `**` for any depth). Quote the pattern when
  it has spaces. With a wildcard, give `--out`: by default the output goes into
  the first matched folder.
- `"D:\Astro\"` with a trailing backslash makes Windows pass a stray `"`
  (and can glue the next options onto the path). The script repairs this, but
  it is cleaner to leave the trailing backslash off.
- `-j` is capped at 61 on Windows (a Python limit). Each worker needs about
  1 GB of RAM for a 26 MP frame.
- Output is plain text. Characters the console cannot show (for example a
  non-ASCII filter name on a cp1252 console) are printed as `?`. The CSV is
  UTF-8 with a BOM so Excel opens it correctly.

Folders are searched recursively. Only frames with `IMAGETYP` = Light are used,
so flats and darks in the same tree are ignored. The output folder is skipped
when scanning, so re-running is safe. Measurements are cached in
`~/.cache/subselect` (Windows: `%LOCALAPPDATA%\subselect`; set
`SUBSELECT_CACHE` to use another folder), so a re-run with different limits
finishes almost instantly. Copies keep the original timestamps, and a file that
already exists at the destination with the same size and time is not copied
again.

## What it measures

| column | meaning |
|---|---|
| FWHM, minor | median star FWHM and minor axis, arcsec |
| ecc | median eccentricity (0 = round). Trailing, tilt and bad guiding raise it |
| trail | smear shared by all stars (tracking/flexure), arcsec. Random seeing cancels out |
| stars | detected stars. Clouds, haze and bad focus lower it |
| bg x | sky level relative to the same filter, exposure and night (moon, cloud, stray light) |
| grad x | large-scale gradient relative to the same filter, exposure and night (stray light) |
| sky e-, eff | with `--darks`: sky above dark in electrons, and the share of the ideal SNR you keep after read noise, `sky/(sky+RN²)` |

## Grades

Each check has a "good" limit and a "reject" limit. Anything in between is
borderline:

| check | good if | reject if | flags |
|---|---|---|---|
| FWHM | ≤ 4.5″ | > 6.0″ | `--fwhm-good`, `--fwhm-max` |
| eccentricity | ≤ 0.72 | > 0.82 | `--ecc-good`, `--ecc-max` |
| stars vs median | ≥ 0.5× | < 0.25× | `--stars-good`, `--stars-min`, `--min-stars` |
| background vs median | ≤ 1.5× | > 2× | `--bg-good`, `--bg-max` |
| gradient vs median | ≤ 2× | > 3× | `--grad-good`, `--grad-max` |

The FWHM limits depend on your seeing and image scale. Run once in report mode
and set them from your own best subs.

## Output layout

```
_selected/NIGHT_2026-10-27/LIGHT/       good subs, all filters in one folder
_selected/NIGHT_2026-10-27/FLAT/        that night's flats, for the filters used
_selected/NIGHT_2026-10-27/DARKFLAT/    that night's flat-darks, if any
_selected_borderline/NIGHT_2026-10-23/LIGHT/
```

If a night and filter has a master flat in the source tree (IMAGETYP
`Master Flat`, e.g. built by [masterflats](../masterflats)), only the master is
copied into `FLAT/`, and no raw flats or flat-darks for it. If an earlier run
already copied raw flats there, the script says so; delete them so WBPP uses
only the master. Any folder whose name starts with `_selected` is never
scanned, so other targets' output trees don't get picked up.

Night is the evening date. It is taken from a `NIGHT_YYYY-MM-DD` folder in the
source path, or from local time minus 12 h, so dawn flats belong to the night
before. At the end, the script lists any night and filter that has lights but no
flats.

In PixInsight WBPP:
1. Add Directory: `_selected`. Also add `_selected_borderline` for a
   comparison run. WBPP sorts lights and flats by type and filter itself.
2. Grouping Keywords: add `NIGHT`, tick **Pre-processing**, untick
   **Post-processing**. Each night's lights then get that night's flats, and
   all nights still stack into one master per filter.

## Tests

```
pip install pytest
python -m pytest subselect -q        # from the repository root
```

The tests build small synthetic XISF and FITS subs in a temporary folder and
run on Linux and Windows.
