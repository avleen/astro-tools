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
| [flextest](flextest) | Diagnoses flexure and tracking: in-sub trails, drift and field rotation between subs, compensator on/off comparisons, sky background vs read noise. |
| [masterflats](masterflats) | Turns each night's raw flats into one master flat per filter, next to the raw flats. Subtracts flat-darks, a master bias or a master dark. Writes zstd-compressed XISF that WBPP uses as-is. Idempotent, so it can run every day; can delete the raw flats afterwards. |
| [nina-sequences](nina-sequences) | N.I.N.A. Advanced Sequencer files: an unattended nightly sequence built on Target Scheduler, ending with the flats that masterflats uses. |
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

To run masterflats and subselect together each morning, see
[Workflow example](#workflow-example).

Quote paths that contain spaces. Wildcards such as `*.xisf` work in cmd and
PowerShell too: the tools expand them themselves.

The output layout is `_selected/NIGHT_<date>/LIGHT` plus `.../FLAT`. Add the
folder in WBPP and add `NIGHT` as a Pre-processing-only grouping keyword. Each
night's lights are then calibrated with that night's flats. See
[subselect/README.md](subselect/README.md) for the limits and all the options.

## Notes

- XISF can be uncompressed or lz4, zstd or zlib compressed, with or without byte shuffling. FITS must be uncompressed.
- Measurements are cached, so re-running with different limits is instant.
- Expect about 1 GB of RAM per worker for a 26 MP frame (`-j` sets the number of workers).
- Tested with an ASI2600MM Pro and N.I.N.A. 3.x files, on Windows and Linux.
- Tests: `pip install pytest`, then `python -m pytest` from the repo root.

## Workflow example

A night's data from the rig, as N.I.N.A. saves it:

```
D:\Astro\Downloads\Heart Nebula\
    NIGHT_2026-10-27\
        LIGHT\H\  LIGHT\O\  LIGHT\S\     300 s subs
        FLAT\H\   FLAT\O\   FLAT\S\      25 flats per filter, taken at dawn
```

**1. Build the master flats.** Run with `-n` first to see what it will do, then
again without it:

```
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master" -n
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master"
```

Each `FLAT\<filter>\` folder now also holds `MasterFlat_2026-10-27_<filter>_<exp>s.xisf`.
Earlier nights whose masters are up to date are skipped, so pointing it at
the whole downloads folder every day is fine.

**2. Grade the subs.** A report only; nothing is copied or moved:

```
python subselect/subselect.py "D:\Astro\Downloads\Heart Nebula" --darks "D:\Astro\Darks -10C\master"
```

Check the rejects. If the limits are too strict or too loose, change them
(`--fwhm-good`, `--ecc-good` and so on) and run it again. Measurements are
cached, so this is instant.

**3. Copy the keepers:**

```
python subselect/subselect.py "D:\Astro\Downloads\Heart Nebula" --copy
```

```
D:\Astro\Downloads\Heart Nebula\
    _selected\
        NIGHT_2026-10-27\
            LIGHT\...                     good subs
            FLAT\...                      the master flats, not the raw flats
    _selected_borderline\
        NIGHT_2026-10-27\LIGHT\...        borderline subs
```

**4. Stack.** In WBPP, add `_selected` (and `_selected_borderline` if you want
them), add `NIGHT` as a Pre-processing-only grouping keyword, and add your
master darks. Each night's lights are calibrated with that night's master
flats, and all nights stack together.

**5. Tidy up (optional).** Once you're happy with the masters, delete the raw
flats:

```
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master" --delete-raw
```

The masters stay where they are and are still picked up by later runs of either tool.

## Logo

In [docs/logo](docs/logo): SVGs to scale to any size, PNGs from 16 to 1024 px,
and a favicon.

## License

MIT
