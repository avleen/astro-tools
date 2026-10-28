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
- **Cache:** star lists are cached in `cache/` next to the script.

## skybg.py

Measures sky background above a master dark, in electrons, and how well it
swamps read noise.

```
python skybg.py MASTER_DARK LIGHT...
```

Efficiency is `sky / (sky + RN²)`, the share of the ideal SNR you keep. Read
noise is fixed at 1.2 and 1.5 e⁻ (ASI2600MM at gain 100). Gain comes from the
`EGAIN` header. For a calibrated sky report with other cameras, use
`subselect.py --darks ... --read-noise N` instead.
