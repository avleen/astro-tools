# astro-tools

Python tools for checking and sorting deep-sky subs from a remote imaging rig.
They read N.I.N.A. XISF directly, including lz4/zstd compressed files, and
uncompressed FITS. No PixInsight is needed to run them.

| tool | what it does |
|---|---|
| [subselect](subselect) | Grades every light sub (FWHM, eccentricity, trailing, star count, background, gradient, optional sky e⁻ and read-noise efficiency). Copies the good and borderline ones, plus the same night's flats, into a layout that drops straight into PixInsight WBPP. |
| [flextest](flextest) | Diagnoses flexure and tracking: in-sub trails, drift and field rotation between subs, compensator on/off comparisons, sky background vs read noise. |

## Install

```
git clone https://github.com/avleen/astro-tools.git
cd astro-tools
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt        # Windows: .venv\Scripts\pip install -r requirements.txt
```

## Typical use

```
# grade a target's subs; report only, nothing is changed
.venv/bin/python subselect/subselect.py "D:/Astro/Heart Nebula" --darks "D:/Astro/Darks -10C"

# happy with the grades: build _selected/ and _selected_borderline/
.venv/bin/python subselect/subselect.py "D:/Astro/Heart Nebula" --copy
```

The output layout is `_selected/NIGHT_<date>/LIGHT` plus `.../FLAT`. Add the
folder in WBPP and add `NIGHT` as a Pre-processing-only grouping keyword. Each
night's lights are then calibrated with that night's flats. See
[subselect/README.md](subselect/README.md) for the limits and all the options.

## Notes

- Measurements are cached, so re-running with different limits is instant.
- Expect about 1 GB of RAM per worker for a 26 MP frame (`-j` sets the number of workers).
- Tested with an ASI2600MM Pro and N.I.N.A. 3.x files.

## License

MIT
