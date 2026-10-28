# subselect

Grades astrophotography light subs (XISF from N.I.N.A., or uncompressed FITS) on
star size, shape, star count and background. It can then copy the keepers into a
folder tree that PixInsight WBPP can pair with same-night flats.

## Install

```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt      # Windows: .venv\Scripts\pip ...
```

## Use

```
python subselect.py "D:/Astro/Heart Nebula"                 # report only, changes nothing
python subselect.py "D:/Astro/Heart Nebula" --copy          # copy to "D:/Astro/Heart Nebula/_selected"
python subselect.py lights/ --darks masters/ --csv out.csv  # add sky e- / read-noise efficiency, CSV
```

Folders are searched recursively. Only frames with `IMAGETYP` = Light are used,
so flats and darks in the same tree are ignored. The output folder is skipped
when scanning, so re-running is safe. Measurements are cached in
`~/.cache/subselect`, so a re-run with different limits finishes almost
instantly.

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
