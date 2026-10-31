<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/logo/logo-mark-dark.svg">
    <img src="docs/logo/logo-mark-light.svg" width="180" alt="astro-tools logo">
  </picture>
</p>

<h1 align="center">astro-tools</h1>

<p align="center"><b>Grade subs · build master flats · chase down flexure</b></p>

A collection of small command-line tools for handling XISF and FITS files in
an astrophotography workflow: grading subs, building master flats, measuring
flexure and recovering lost work. Each one does a single job, and they're
built to be chained together.

| tool | what it does |
|---|---|
| [subselect](subselect) | Grades every light sub (FWHM, eccentricity, trailing, star count, background, gradient, optional sky e⁻ and read-noise efficiency). Copies the good and borderline ones, plus the same night's flats, into a layout that drops straight into PixInsight WBPP. |
| [subdepth](subdepth) | Measures what each sub is worth on a faint target: plate-solves with ASTAP, gets a zero point from the Gaia star database, and reports sky brightness (mag/arcsec²), transparency and effective time, per sub, night and filter. |
| [flextest](flextest) | Diagnoses flexure and tracking: in-sub trails, drift and field rotation between subs, compensator on/off comparisons, sky background vs read noise. |
| [masterflats](masterflats) | Turns each night's raw flats into one master flat per filter, next to the raw flats. Subtracts flat-darks, a master bias or a master dark. Writes zstd-compressed XISF that WBPP uses as-is. Idempotent, so it can run every day; can delete the raw flats afterwards. |
| [nina-sequences](nina-sequences) | N.I.N.A. Advanced Sequencer files: an unattended nightly sequence built on Target Scheduler, ending with the flats that masterflats uses. |
| [processing](processing) | PixInsight notes: from WBPP masters to a finished narrowband nebula with RGB stars, with settings for the steps that need care (MGC, NarrowbandNormalization, LHE, star reduction). |
| [pi-swap-recover](pi-swap-recover) | Recovers images from PixInsight's swap files after a crash, converting every undo-history state to XISF, with optional previews. |

## Install

Python 3.9 or newer, on Windows, Linux or macOS.

Windows (PowerShell or cmd), using the `py` launcher from python.org:

```
git clone https://github.com/avleen/astro-tools.git
cd astro-tools
py -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

Linux / macOS / WSL:

```
git clone https://github.com/avleen/astro-tools.git
cd astro-tools
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

`pillow` is only needed for `pi-swap-recover --previews`; install it with `pip install pillow`.

subdepth also needs [ASTAP](https://www.hnsky.org/astap.htm) with its H18 or
H17 star database.

## Quick start

With the venv activated, after each night:

```
python masterflats/masterflats.py "D:\Astro" --darks "D:\Astro\Darks"
python subselect/subselect.py "D:\Astro\Heart Nebula" --copy
python subdepth/subdepth.py "D:\Astro\Heart Nebula\_selected"
```

When `eff h` per filter reaches your goal, stack `_selected` in WBPP and process.

## Workflow

### Folder layout

The tools work on folders of XISF or FITS files and read everything they need
from the headers (`IMAGETYP`, `FILTER`, `EXPTIME`, temperature, gain). A layout
like N.I.N.A.'s works well:

```
D:\Astro\Heart Nebula\
    NIGHT_2026-10-27\
        LIGHT\H\  LIGHT\O\  LIGHT\S\
        FLAT\H\   FLAT\O\   FLAT\S\
```

The night is taken from a `NIGHT_<date>` folder in the path, or else from the
local time minus 12 h. Folders are searched recursively. Quote paths with
spaces. Wildcards such as `*.xisf` work in cmd and PowerShell too.

### 1. Master flats

```
python masterflats/masterflats.py "D:\Astro" --darks "D:\Astro\Darks" -n
python masterflats/masterflats.py "D:\Astro" --darks "D:\Astro\Darks"
```

`-n` shows what would be built. Each night and filter gets one master flat,
next to its raw flats. Nights whose masters are up to date are skipped, so it
can run on the whole data folder every day.

| change | option |
|---|---|
| flat-darks, a master bias, or master darks for the flats | `--darks` (repeatable; matched by exposure and temperature) |
| nights with few flats | `--min-flats` (default 3) |
| delete raw flats once their master exists | `--delete-raw` |

### 2. Grade the subs

```
python subselect/subselect.py "D:\Astro\Heart Nebula"
```

Report only. Check the rejects and borderline subs, adjust the limits if
needed, and run again. Measurements are cached.

| change | option |
|---|---|
| FWHM limits, arcsec | `--fwhm-good`, `--fwhm-max` (4.5, 6.0) |
| elongation | `--ecc-good`, `--ecc-max` (0.72, 0.82) |
| star count vs median | `--stars-good`, `--stars-min` (0.5, 0.25) |
| background, gradient vs the night's median | `--bg-good`, `--bg-max`, `--grad-good`, `--grad-max` |
| also report sky e⁻ and read-noise efficiency | `--darks`, `--read-noise` |
| headers without EGAIN or pixel scale | `--gain`, `--scale` |

### 3. Copy the keepers

```
python subselect/subselect.py "D:\Astro\Heart Nebula" --copy
```

```
D:\Astro\Heart Nebula\
    _selected\NIGHT_2026-10-27\LIGHT\...      good subs
    _selected\NIGHT_2026-10-27\FLAT\...       that night's master flats
    _selected_borderline\NIGHT_.../LIGHT\...  borderline subs
```

`--good-only` skips the borderline subs. `--out` sets another destination.
A later run doesn't remove files, so after changing limits a sub can be in
both folders.

### 4. Measure the depth

```
python subdepth/subdepth.py "D:\Astro\Heart Nebula\_selected"
```

Prints, per filter and night and per filter: hours, effective hours (`eff h`),
sky brightness (`mu`), transparency (`T`) and zero point. `eff h` per filter is
the progress measure. Low `T` points to cloud or haze, low `eff %` to a bright
sky (moon, twilight).

Each sub is plate-solved once and cached. The reference "ideal sub" is the best
of the subs in the run, so earlier nights' figures can shift slightly as new
nights are added.

| change | option |
|---|---|
| read noise of your camera at its gain, e⁻ | `--read-noise` (default 1.4) |
| bias level for a camera without an `OFFSET` estimate | `--darks` or `--pedestal` |
| frames without EGAIN | `--gain`, or `--gain 100=0.25 --gain 0=0.79` per gain setting |
| ASTAP in another folder | `--astap` |
| per-sub table | `-v`, `--csv` |

Check the `sky from` line and any bias/noise disagreement report: they mean
the gain, read noise or bias level is wrong for those nights. See
[subdepth/HOW-IT-WORKS.md](subdepth/HOW-IT-WORKS.md).

### 5. Decide whether there is enough

Depth grows with the square root of effective time: twice the `eff h` reaches
0.38 mag fainter, four times reaches 0.75 mag fainter. Set a goal in `eff h`
per filter, and stop when the next doubling isn't worth the nights it would
take.

### 6. Stack

In WBPP, add `_selected` (and `_selected_borderline` if wanted), add `NIGHT` as
a Pre-processing-only grouping keyword, and add your master darks. Each
night's lights are calibrated with that night's master flats, and all nights
stack together. PSF Signal Weight weights subs much as effective time does.

### 7. Process

See [processing](processing) for a narrowband workflow from WBPP masters to a
finished image.

### 8. Tidy up (optional)

Once the masters are confirmed good, delete the raw flats:

```
python masterflats/masterflats.py "D:\Astro" --darks "D:\Astro\Darks" --delete-raw
```

Later runs of either tool still find the masters.

## Notes

- XISF can be uncompressed or lz4, zstd or zlib compressed, with or without byte shuffling. FITS must be uncompressed.
- Measurements are cached; re-running with different limits doesn't re-measure.
- Expect about 1 GB of RAM per worker for a 26 MP frame (`-j` sets the number of workers).
- Tested with an ASI2600MM Pro and N.I.N.A. 3.x files, on Windows and Linux.
- Tests: `pip install pytest`, then `python -m pytest` from the repo root.

## Logo

In [docs/logo](docs/logo): SVGs to scale to any size, PNGs from 16 to 1024 px,
and a favicon.

## License

MIT
