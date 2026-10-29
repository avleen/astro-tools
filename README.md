# astro-tools

Python tools for checking and sorting deep-sky subs from a remote imaging rig.
They read N.I.N.A. XISF directly, including lz4/zstd compressed files, and
uncompressed FITS. No PixInsight is needed to run them.

| tool | what it does |
|---|---|
| [subselect](subselect) | Grades every light sub (FWHM, eccentricity, trailing, star count, background, gradient, optional sky e⁻ and read-noise efficiency). Copies the good and borderline ones, plus the same night's flats, into a layout that drops straight into PixInsight WBPP. |
| [flextest](flextest) | Diagnoses flexure and tracking: in-sub trails, drift and field rotation between subs, compensator on/off comparisons, sky background vs read noise. |
| [masterflats](masterflats) | Turns each night's raw flats into one master flat per filter, next to the raw flats. Subtracts flat-darks, a master bias or a master dark. Idempotent, so it can run every day; can delete the raw flats afterwards. |
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

## Typical use

With the venv activated (`python` below is the venv's interpreter on every platform):

```
# grade a target's subs; report only, nothing is changed
python subselect/subselect.py "D:\Astro\Heart Nebula" --darks "D:\Astro\Darks -10C"

# happy with the grades: build _selected/ and _selected_borderline/
python subselect/subselect.py "D:\Astro\Heart Nebula" --copy
```

Each tool does one job, and they work together. A daily routine:

```
# 1. new flats -> master flats (already-built masters are skipped)
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master"

# 2. grade and copy; subselect copies the masters instead of the raw flats
python subselect/subselect.py "D:\Astro\Downloads\Heart Nebula" --copy
```

Quote paths that contain spaces. Wildcards such as `*.xisf` work in cmd and
PowerShell too: the tools expand them themselves.

The output layout is `_selected/NIGHT_<date>/LIGHT` plus `.../FLAT`. Add the
folder in WBPP and add `NIGHT` as a Pre-processing-only grouping keyword. Each
night's lights are then calibrated with that night's flats. See
[subselect/README.md](subselect/README.md) for the limits and all the options.

## Notes

- Measurements are cached, so re-running with different limits is instant.
- Expect about 1 GB of RAM per worker for a 26 MP frame (`-j` sets the number of workers).
- Tested with an ASI2600MM Pro and N.I.N.A. 3.x files, on Windows and Linux.
- Tests: `pip install pytest`, then `python -m pytest` from the repo root.

## License

MIT
