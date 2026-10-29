<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# flextest

These scripts diagnose tracking and flexure problems. They don't pick subs;
use [subselect](../subselect) for that.

## flextest.py

Measures star trailing inside each sub and image drift between subs. When the
two agree, the main camera is moving relative to the guide scope: differential
flexure.

```
python flextest.py FILE...                     # one sequence
python flextest.py --map FILE...               # 3x3 field maps of FWHM/ecc/trail direction
python flextest.py --off F1 F2... --on F3 F4...  # e.g. flexure compensator off vs on
```

- **Dithering:** turn it off for the test subs, or the between-sub drift numbers mean nothing.
- **How many subs:** six 60 s subs per block measure drift to about ±0.04″/min.
- **Field rotation:** the "rot" column shows it. Rotation that grows through the night points to polar alignment error.
- **Cache:** star lists are cached in `cache/` next to the script, keyed on the
  file's header and size, so the same sub hits the cache from Windows, WSL or a
  copy on another drive. Set `FLEXTEST_CACHE` to use another directory. If the
  directory can't be written, the script warns and carries on uncached.

## skybg.py

Measures sky background above a master dark, in electrons, and how well it
swamps read noise.

```
python skybg.py MASTER_DARK LIGHT...
```

Both scripts expand wildcards themselves, so `*.xisf` works in cmd.exe and
PowerShell as well as in bash.

Efficiency is `sky / (sky + RN²)`, the share of the ideal SNR you keep. Read
noise is fixed at 1.2 and 1.5 e⁻ (ASI2600MM at gain 100). Gain comes from the
`EGAIN` header. For a calibrated sky report with other cameras, use
`subselect.py --darks ... --read-noise N` instead.

## Windows

Native Windows (cmd.exe or PowerShell) works as well as Linux and WSL. Set up a
venv once from the repo root:

```
py -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

In PowerShell, activate with `.venv\Scripts\Activate.ps1`.

Then run the scripts from the repo root:

```
python flextest\flextest.py "D:\Astro\Heart Nebula\NIGHT_2026-10-27\LIGHT\*_H_*.xisf"
python flextest\flextest.py --off "D:\data\off\*.xisf" --on "D:\data\on\*.xisf"
python flextest\skybg.py "D:\Astro\Darks\master_dark_600s.xisf" "D:\Astro\Heart Nebula\NIGHT_2026-10-27\LIGHT\*.xisf"
```

- **Quoting:** put paths that contain spaces in double quotes. A wildcard
  inside quotes still works, because the scripts expand it themselves, not
  the shell.
- **Wildcards:** `*`, `?` and `[...]` work in any path component, for example
  `D:\data\*\LIGHT\H\*.xisf`. Use `**` to match any depth of folders. A
  pattern that matches nothing is an error.
- **Console output:** characters the console can't show print as `?` rather
  than crashing. This also applies when output is redirected to a file.
- **Cache directory:** `set FLEXTEST_CACHE=C:\Temp\flexcache` in cmd.exe, or
  `$env:FLEXTEST_CACHE = "C:\Temp\flexcache"` in PowerShell.

## Tests

```
pip install pytest
python -m pytest flextest -q
```

The tests build small synthetic XISF subs: drifting, trailed stars and round
static stars, uncompressed and lz4/zstd compressed. They check XISF reading,
trail measurement, drift recovery, wildcard expansion and redirected cp1252
output. They run on Linux and Windows.
