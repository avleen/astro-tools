<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# subdepth

Measures sky brightness, transparency and effective time for each light sub,
using ASTAP plate solving and the Gaia star catalogue.

Effective time is the exposure of an ideal sub (transparency 1, the darkest sky
of that filter in the run, no read noise) with the same signal-to-noise on a
target fainter than the sky.

Method, formulas and measured constants: [HOW-IT-WORKS.md](HOW-IT-WORKS.md).

## Install

Same as [subselect](../subselect) (numpy, scipy, lz4, zstandard), plus
[ASTAP](https://www.hnsky.org/astap.htm) with the **H18** or **H17** star
database. Those databases are an extract of Gaia (BP magnitudes). The Gaia data
is free to use with credit to ESA/Gaia/DPAC.

subdepth looks for ASTAP in `$ASTAP_DIR`, then `C:\Program Files\astap`
(also from WSL), `/opt/astap` and `/usr/share/astap`. Use `--astap` for any
other folder.

## Use

```
python subdepth.py "D:\Astro\Heart Nebula" --darks "D:\Astro\Calibration Frames\Askar FRA600\Darks - -10C, Gain 100\master"
python subdepth.py "D:\Astro\Heart Nebula" --read-noise 1.25 --csv depth.csv -v
python subdepth.py "D:\Astro\Dark Shark\_selected" --gain 100=0.25 --gain 0=0.79
python subdepth.py "D:\Astro\Old Target" --pedestal 2750 --gain 0.25
```

- Folders are searched recursively; only `IMAGETYP` = Light frames are used.
- ASTAP solves a temporary FITS copy of each sub in the temp folder. Under WSL
  that is the Windows temp folder, because Windows ASTAP can't open
  `\\wsl.localhost` paths.
- Measurements are cached in `~/.cache/subdepth` (Windows:
  `%LOCALAPPDATA%\subdepth`; set `SUBDEPTH_CACHE` to move it). A re-run with
  different darks, pedestal, gain or read noise doesn't re-solve.
- One-shot colour frames (with a `BAYERPAT` header) are measured as 2×2
  superpixels: R+G+G+B summed into one mono pixel.
- `-v` prints one line per sub; `--csv` writes them all.

## What it measures

| column | meaning |
|---|---|
| zp | zero point: the Gaia BP magnitude of a star that gives 1 e⁻/s. From stars of magnitude 10–15 matched to the catalogue |
| sky e- | sky above the bias level, per pixel per sub, in electrons |
| b/n | sky from the bias level ÷ sky from the pixel noise. Close to 1 when gain, read noise and bias are right |
| gain? | the gain implied by the sky level and its noise. Should be close to EGAIN |
| mu | sky brightness, mag/arcsec², on the Gaia BP scale |
| T | transparency: star flux compared with the clearest subs of that filter (haze, thin cloud, dew, low altitude) |
| eff | read-noise efficiency, `sky / (sky + RN²)` |
| t_eff | effective time in seconds (below) |

**Sky.** Measured on the darkest 20% of 128-px tiles in the central 60% of the
frame, with stars masked, two ways:

- From the bias level: sub − bias. The bias comes from a matching master dark
  (`--darks`), else `--pedestal`, else an estimate from the `OFFSET` header
  (ASI2600MM and ASI533MC only, for now).
- From the pixel noise: noise² = sky + read noise². This needs no bias level.

`b/n` is the ratio of the two. When a night's median ratio is outside
1/1.1–1.1, subdepth reports it: the gain, read noise, offset or darks are
probably wrong for those subs. An estimated bias that fails this check is not
used for that night; the noise value is used instead.

**Darks** are optional. Only a master dark's median (its bias level) is used.
See [Darks and the sky measurement](HOW-IT-WORKS.md#darks-and-the-sky-measurement).

**Gain.** `EGAIN` ÷ the storage step. 14-bit cameras such as the ASI533MC store
values ×4 in 16-bit files, while EGAIN is per native ADU. subdepth detects the
step and says when it divides. Frames without EGAIN need `--gain`, in e⁻ per
file ADU: one value, or one per camera gain setting
(`--gain 100=0.25 --gain 0=0.79`, matched on the `GAIN` header).

**Read noise.** `--read-noise` is per sub in e⁻ (default 1.4). It changes with
camera gain: 1.2–1.3 e⁻ on the ASI2600MM at gain 100, about 3.5 e⁻ on the
ASI533MC at gain 0.

**Effective time.**

```
t_eff = t × T² × b_ref / (b + RN²/t)        b = sky rate, e-/s/px
```

The references (95th percentile of zp, 5th percentile of sky rate) are taken
per filter from the subs in the run, so compare datasets within one run.
t_eff is capped at the sub's own exposure.

**mu** is on the Gaia BP scale for every filter. It is comparable night to
night within a filter and camera. For narrowband it is not comparable to an
SQM reading.

## Output

```
filter night              subs  hours  eff h eff %     mu     T  sky e-     zp
H 2026-10-25                 3   0.50   0.06    13  18.00  0.99   102.3  15.47
H 2026-10-27                10   1.67   0.44    27  18.36  0.97    65.8  15.40
H 2026-10-28                 8   0.67   0.03     5  17.13  0.30    37.1  14.13
```

`eff h` is the sum of t_eff in hours. Here 28 Oct lost 70% of its starlight
(T = 0.30), so 40 minutes of Ha were worth 2 minutes.

Subs that ASTAP can't solve (badly trailed or clouded out) get no zero point
and no t_eff; they're listed at the top.
