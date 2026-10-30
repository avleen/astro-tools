# astro-tools

Small command-line tools for handling XISF and FITS files in an
astrophotography workflow. Each tool is one script in its own folder with its
own README and tests.

| folder | script | job |
|---|---|---|
| `subselect/` | `subselect.py` | grade light subs, copy keepers + that night's flats into `_selected/NIGHT_<date>/{LIGHT,FLAT}` |
| `masterflats/` | `masterflats.py` | one master flat per night and filter, next to the raw flats; idempotent |
| `flextest/` | `flextest.py`, `skybg.py` | flexure, tracking and sky-background diagnostics |
| `pi-swap-recover/` | `pi_swap_recover.py` | recover images from PixInsight swap files |
| `nina-sequences/` | `*.json` | N.I.N.A. Advanced Sequencer files (data, not code; no tests) |
| `processing/` | `README.md` | PixInsight processing notes (docs only; no tests) |
| `docs/logo/` | `make_logo.py` | generates every logo file in that folder |

## Design rules

- Lots of small tools working together, not one big app. A new job gets a new
  script, not a new mode on an existing one.
- Dependencies: numpy, scipy, lz4, zstandard (`requirements.txt`). pillow is
  optional (pi-swap-recover previews only). Don't add more without asking.
- `masterflats` imports its file readers from `../subselect/subselect.py` via
  `sys.path`; keep shared readers there rather than copying them.
- Everything must work on native Windows as well as Linux/WSL: quote-safe
  paths, wildcard expansion done by the tools themselves (`expand`, `fix_argv`),
  no shell-only tricks.
- Never change or delete a user's raw subs or flats unless a flag explicitly
  says so (`--copy`, `--delete-raw`). Masters and caches are written to a
  `.tmp` file and `os.replace`d. Report-only (or `-n`) is the default or
  always available.
- XISF output is Float32, zstd + byte shuffling (`compression="zstd+sh:..."`),
  so PixInsight and WBPP open it as-is.
- Night = the `NIGHT_YYYY-MM-DD` folder in the path, else local time minus 12 h.

## Tests

From the repo root, with both venvs:

```
./.venv/bin/python -m pytest -q                  # Linux / WSL
./.venv-win/Scripts/python.exe -m pytest -q      # native Windows
```

Both must pass before committing. New behaviour gets a test in that tool's
`tests/` folder; tests build small synthetic XISF/FITS files, never real data.

## Docs

- Each tool's README, the root README and the script's `--help` docstring
  must agree. Update all three when options or behaviour change.
- Write plainly: short sentences, concrete examples, real paths. No marketing
  language or filler.
- Examples use Windows paths (`D:\Astro\...`), since that's where the rig's
  data lives.

## Git

Commit or push only when asked. Work on main unless asked to branch.
