# Depth measurement: design

2026-10-31. subselect grades sharpness. This adds depth: how much signal-to-noise
each sub, and each stack, gives on a faint target.

## Tools

| tool | input | output | status |
|---|---|---|---|
| `subdepth` | raw light subs (+ master darks) | per sub: zero point, sky e⁻ and mag/arcsec², transparency, effective time | built |
| `stackdepth` | a master (+ its registered `_r` subs) | surface brightness limit; noise vs N; odd/even split noise | next |
| `sitecompare` | subdepth CSVs or SQM values for two sites | hours here = hours there, per filter | later |

Python only. A PixInsight script may come later for marking regions (the one
interactive step); measuring stays on the command line.

## Catalogue

ASTAP's local H18/H17 star database (Gaia eDR3, BP magnitude, to mag 17/18).
No remote calls. The `.1476` format: 36 declination rings, each split into
equal RA cells; 110-byte header (byte 109 = record size 5); 5-byte records
(3-byte RA, low 2 bytes of Dec); a record with RA = 0xFFFFFF starts a group and
carries magnitude×10 + 16 and the Dec high byte + 128. Read from ASTAP's
source (`unit_command_line_star_database.pas`) and checked against the files.

ASTAP CLI 2022 reads FITS only (not XISF) and has no `-sqm`, so subdepth writes
a temporary float32 FITS and solves that. Windows ASTAP under WSL can't open
`\\wsl.localhost` paths, so the temp file goes in the Windows `%TEMP%`.

## Formulas

- Zero point: median over matched mag 10–15 stars of `mag + 2.5 log10(flux e⁻/s)`.
  Aperture 1.5 × FWHM, annulus background.
- Sky: darkest 20% of 128-px tiles in the central 60%, stars masked.
  - from the bias: `(sky − bias) × gain`; bias = matching master dark, else
    `--pedestal`, else `k × OFFSET × storage step` (k = 10.0 ASI2600MM, 9.8
    ASI533MC native ADU per offset unit, from the dark comparison)
  - from the noise: `(σ × gain)² − RN²`, σ from neighbouring-pixel differences
  - always both when a bias is known: `b/n` ratio and `gain? = s/σ²` (solved
    with RN). An estimated bias is used per night only if the median b/n is
    within 10%; otherwise that night uses the noise value. Disagreements are
    reported.
- `mu = zp − 2.5 log10(sky e⁻/s/arcsec²)`
- `T = 10^(−0.4 (zp_ref − zp))`, zp_ref = 95th percentile per filter
- `t_eff = t · T² · b_ref / (b + RN²/t)`, b_ref = 5th percentile sky rate per filter
- stackdepth: `μ_lim = zp − 2.5 log10(3 σ_box / scale²)`, σ_box = scatter of
  10″ box means on empty sky (measured, so correlated noise counts)

## Validation so far (2026-10-31)

- Sky by dark and by noise agree within 5% on 2600MM subs; `gain?` 0.24–0.25 vs
  EGAIN 0.243. B 60 s sky 17.7 e⁻ vs 17.1 e⁻ measured earlier with subselect.
- Moonlit vs dark Ha 13% and SII 9% effective, matching the earlier estimate of
  1/7 and 1/12.
- ASI533MC: N.I.N.A. writes EGAIN per native 14-bit ADU (1.0056 at gain 100,
  3.18 at gain 0) but stores values ×4. Using EGAIN as-is made the N.I.N.A.
  nights' zero point 1.4 mag off the non-N.I.N.A. night (whose ADU zero point
  agreed). Fix: detect the storage step from the pixel values, divide EGAIN.
  File-scale gain is then 0.25 (gain 100) and 0.79 (gain 0). An estimate of the
  zero point from aperture and bandwidth was 1.8 mag off, so cross-night
  consistency is the better check than theory.

- Bias from OFFSET vs noise (no darks), 2026-10-31: 2600MM all six filters
  b/n 1.01–1.05; 533MC N.I.N.A. gain 100 1.04; 533MC gain 0 about 0.9 (EGAIN
  3.18 or the bias slightly off); frames from a non-N.I.N.A. capture program
  are flagged (different bias, no storage step).

## Open

- Absolute mu not yet checked against an independent SQM measurement (ASTAP
  GUI: Tools → SQM, on the same sub).
- Badly trailed subs don't solve, so they get no zero point.
- stackdepth and sitecompare not started.
