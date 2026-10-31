<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# How subdepth works

This describes each step, the formulas, and the measurements behind the
defaults. For options and examples, see the [README](README.md).

## Steps

For each light sub:

1. Read the frame. Detect the storage step. Sum colour frames into 2×2 superpixels.
2. Measure the sky level and the per-pixel noise in star-free sky.
3. Find stars, plate-solve with ASTAP, match stars to the Gaia catalogue, and
   measure the zero point.
4. Cache those measurements, which are all in ADU.

Then, for the whole run:

5. Convert to electrons using the gain, the read noise and the bias level.
6. Compute sky brightness, transparency and effective time.
7. Summarise per filter and night, and per filter.

Steps 1–4 are cached. Steps 5–7 run every time, so a change to `--gain`,
`--read-noise`, `--darks` or `--pedestal` uses the cached measurements.

## 1. Reading the frame

- **Storage step.** 14-bit cameras write their values ×4 into 16-bit files, so
  every value is a multiple of 4. subdepth samples about a million pixels. If
  every non-zero integer value is a multiple of 16, 8, 4 or 2, that is the
  storage step. Float data has step 1.
- **Colour frames** (a `BAYERPAT` header) are summed in 2×2 blocks (R+G+G+B)
  into one mono pixel. In electrons that makes:
  - pixel scale ×2
  - read noise ×2 (four independent read noises: √4)
  - bias level ×4

## 2. Sky level and noise

The sky is measured on the central 60% of the frame, which avoids corner
vignetting and amp glow.

1. Split it into 128 × 128 px tiles and subtract each tile's median.
2. Mask stars: pixels more than 2σ above the background after a 3 × 3 smoothing,
   grown by 4 px.
3. For each tile, take:
   - the median of its unmasked pixels (the level)
   - σ from differences between horizontal neighbours: 1.4826 × MAD(Δ) / √2.
     Neighbour differences cancel gradients and any nebula structure larger than
     a pixel, so σ is pixel noise only.
4. Keep the darkest 20% of tiles by level, which keeps out most nebula.
5. The sky level and σ are the medians over those tiles.

In a frame filled with nebula, the darkest tiles still contain some of it, so
the sky includes it.

## 3. Zero point

**Plate solve.** The (superpixel) image is written as a temporary float32 FITS
and solved with `astap_cli`:

```
astap_cli -f <temp.fits> -z 0 -fov <height in deg> -ra <hours> -spd <dec+90> -r 10
```

- The pointing comes from `RA`/`DEC`, else from `OBJCTRA`/`OBJCTDEC`. Without
  either, it is a blind solve (`-r 180`).
- The field of view comes from `XPIXSZ` and `FOCALLEN`.
- The WCS is read from ASTAP's `.ini` (`CRPIX`, `CRVAL`, `CD`), and all temporary
  files are deleted.
- Under WSL, the temporary file goes in the Windows `%TEMP%`, because Windows
  ASTAP can't open `\\wsl.localhost` paths.
- Only the TAN projection is used, without distortion terms. That is accurate
  enough to match stars within a FWHM.

**Catalogue.** ASTAP's H18 or H17 database: Gaia eDR3 with BP magnitudes, read
directly from its `.1476` files (format notes are in the script header). Only
the cells covering the field are read.

**Matching and photometry.**

1. Stars are found with subselect's star finder. For colour frames it runs on
   the superpixel image ÷ 4, so its saturation limit applies per raw pixel.
2. Catalogue stars of magnitude 10–15 are projected onto the image. Each
   detection is matched to the nearest one within max(2 px, FWHM). A catalogue
   star claimed by two detections (a close pair) is dropped.
3. Aperture flux: radius max(3 px, 1.5 × FWHM), with the background taken as
   the median of an annulus from r + 4 to r + 12 px. 1.5 × FWHM holds about
   99.8% of a Gaussian star.
4. Per star: `zp = mag + 2.5 log10(flux / t)`. A 3σ clip (σ at least 0.02 mag)
   follows. The zero point is the median, and the scatter is reported as `zp_mad`.
5. At least 10 matched stars are needed. Subs that don't solve or have too few
   matches get no zero point, so no T and no t_eff.

The zero point in ADU becomes one in electrons with `zp = zp_ADU + 2.5 log10(gain)`:
the BP magnitude of a star that gives 1 e⁻/s.

## 4. Electrons

### Gain

`EGAIN` ÷ storage step. Camera drivers write EGAIN per native ADU, so a 14-bit
camera stored ×4 has a quarter of that gain per file ADU. Frames without EGAIN
use `--gain`: either one value, or `SETTING=VALUE` matched on the `GAIN` header.

### Sky in electrons: two methods

**From the bias level:** `sky_bias = (level − bias) × gain`

The bias comes from, in order:
1. a matching master dark (`--darks`: same exposure, sensor temperature within
   2 °C; its median is used)
2. `--pedestal`
3. an estimate from the `OFFSET` header (below)

**From the noise:** `sky_noise = (σ × gain)² − RN²`

This is Poisson statistics: sky variance in electrons equals sky electrons,
plus read noise. It needs no bias level.

The two methods depend on different inputs:

| method | depends on |
|---|---|
| bias | bias level; gain (linear) |
| noise | gain (squared); read noise |

Their ratio (`b/n`) is a check on all of these. With a bias known, the gain
can also be solved from level and noise (the `gain?` column):

```
s = level − bias (ADU),  v = σ² (ADU²)
v·g² = s·g + RN²   →   g = (s + √(s² + 4·v·RN²)) / (2v)
```

### Which sky value is used

- **A master dark or `--pedestal`:** always used. If the night's median b/n is
  outside 1/1.1–1.1, it is reported.
- **An estimated bias:** used for a night only if that night's median b/n is
  within 1/1.1–1.1. Otherwise the whole night uses the noise value, and the
  night is reported. The decision is per night and per estimated value, so a
  night never mixes the two methods.
- **No bias at all:** the noise value.

### Bias from OFFSET

On the cameras measured so far, the bias level is proportional to the
`OFFSET` setting:

```
bias (file ADU per raw pixel) = k × OFFSET × storage step
```

| camera | k (native ADU per offset unit) | measured |
|---|---|---|
| ZWO ASI2600MM | 10.0 | 500.0 ADU at offset 50 |
| ZWO ASI533MC | 9.8 | 196 file ADU at offset 5 (gain 0, darks); about 2,750 at offset 70 (gain 100, inferred from lights) |

The table is `BIAS_PER_OFFSET` in the script. Other cameras get no estimate.
For them, use darks or `--pedestal`, or add the camera once its darks are
measured. The bias also depends on the capture program: frames from a
program that sets a different bias will fail the b/n check.

Checked against the noise method on real subs, with no darks:

| data | b/n |
|---|---|
| ASI2600MM, offset 50, L R G B S H O, 60–1200 s | 1.01–1.05 |
| ASI533MC, gain 100, offset 70, 120 s | 1.04 |
| ASI533MC, gain 0, offset 5, 300 s | about 0.9 |
| ASI533MC, frames from a different capture program | 4.6 (flagged, noise used) |

### Darks and the sky measurement

Master darks of one ASI2600MM and one ASI533MC were compared, in electrons:

| | ASI2600MM −10 °C | ASI2600MM 0 °C | ASI533MC 0 °C, gain 0 |
|---|---|---|---|
| Structure beyond random noise | none measurable | none measurable | none measurable |
| 128-px tile medians, range | 0.08 e⁻ | 0.09 e⁻ | none measurable |
| Corner minus centre (amp glow) | 0.01 e⁻ | 0.01 e⁻ | none measurable |
| Dark current, typical pixel | 0.00016 e⁻/s | 0.0009 e⁻/s | 0.0009 e⁻/s |
| Dark current, mean incl. hot pixels | 0.00035 e⁻/s | 0.0014 e⁻/s | 0.0022 e⁻/s |
| Pixels over 20 e⁻ | 0.18% (1200 s) | 0.11% (180 s) | 0.17% (300 s) |

Apart from hot pixels, these masters are flat at the bias level plus dark
current (ASI2600MM at −10 °C: under 0.5 e⁻ in 20 minutes). The sky is a median
over masked tiles, so hot pixels don't enter it. subdepth uses only a master
dark's median.

Hot pixels in the same darks (not used by subdepth):
- They are linear in exposure: two exposures predict a third within 2–4% per pixel.
- They scale with one factor for temperature: ×2.2–2.3 from −10 °C to 0 °C.
- They are stable: 93–98% of the hottest 0.1% were the same pixels months or
  years apart.

One camera of each model was measured.

## 5. Sky brightness, transparency, effective time

**Sky brightness**

```
rate = sky_e / t / scale²          (e⁻/s/arcsec²)
mu   = zp − 2.5 log10(rate)        (mag/arcsec²)
```

`mu` is on the Gaia BP scale for every filter. It is comparable from night to
night within the same filter and camera. For narrowband filters it is not
comparable to an SQM reading.

**Transparency**

```
T = 10^(−0.4 × (zp_ref − zp)),  capped at 1
```

`zp_ref` is the 95th percentile of `zp` for that filter in the run. T includes
everything that reduces starlight: altitude, haze, thin cloud, dew.

**Effective time.** For a target much fainter than the sky, with target flux
F (e⁻/s per pixel at T = 1), sky rate b (e⁻/s per pixel) and read noise RN:

```
SNR² = (F·T·t)² / (b·t + RN²)
```

An ideal sub has T = 1, sky b_ref and no read noise, so its SNR² = F²·t / b_ref.
Setting the two equal gives the seconds of ideal sub that this sub is worth:

```
t_eff = t × T² × b_ref / (b + RN²/t)
```

- `b_ref` is the 5th percentile of the sky rate for that filter in the run.
- b is floored at b_ref, so t_eff ≤ t.
- SNR² adds when subs are stacked, so t_eff adds too. `eff h` in the summary is
  Σ t_eff / 3600.
- Dark current is left out. On the ASI2600MM at −10 °C it is about 0.2 e⁻ in
  600 s, under 1% of a typical narrowband sky.

Both references come from the subs in the run. t_eff therefore compares subs
within one run, not across separate runs.

`eff` in the per-sub table is the read-noise part alone: `sky / (sky + RN²)`.

## 6. Cache

Per-sub measurements (steps 1–4, all in ADU) are stored as JSON in
`~/.cache/subdepth` (Windows `%LOCALAPPDATA%\subdepth`, or `SUBDEPTH_CACHE`).
They are keyed on the file's path, size and modification time, and the
measurement version. A changed file, or a new version of subdepth, is measured
again.

## Limitations

- Badly trailed or clouded-out subs often don't solve, so they get no zero point.
- Distortion is ignored. With strong distortion, corner stars may not match;
  the zero point then rests on fewer stars.
- The sky in a frame full of nebula includes some nebula.
- mu has not been checked against an independent SQM measurement. Night-to-night
  comparisons don't depend on that.
- `BIAS_PER_OFFSET` is measured on one camera of each model.
