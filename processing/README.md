<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# processing

Notes on processing a narrowband nebula with RGB stars in PixInsight, from
first master to finished image. Written while processing the Heart Nebula
(IC 1805), October 2026, and corrected where steps didn't work as first
described.

Rig: Askar FRA600 (600 mm), ZWO ASI2600MM (IMX571), about 1.29″/px, frame
about 2.2° × 1.5°. Filters H, S, O, R, G, B. Data from
[subselect](../subselect) and [masterflats](../masterflats), stacked in WBPP.

## Overview

The narrowband masters become the nebula and the RGB masters become the stars.
The two are processed separately and only merged at the end.

```
          H   S   O                         R   G   B
          │   │   │                         │   │   │
   ┌──────┴───┴───┴── shared steps, each master (linear) ┴───┴───┴──────┐
   │ merge exposures → DynamicCrop → BXT correct-only → ImageSolver     │
   │ → SPFC → MGC                                                       │
   └──────┬───┬───┬──────────────────────────┬───┬───┬─────────────────┘
          │   │   │                          └─ ChannelCombination → RGB
   NEBULA (starless)                        STARS
   BXT, NXT, StarX → stretch each           SPCC → BXT → StarX (stars)
   → combine → NarrowbandNormalization      → gentle stretch
   → LHE, curves                                         │
          │                                              │
          └──────────── screen stars onto nebula ────────┘
```

Everything stays linear (unstretched) until a step says to stretch. Use an
STF auto-stretch to look at linear images; it only changes the display.

## 1. WBPP

- Master darks at the lights' temperature (−10 °C), and a master flat per
  night (see the [workflow example](../README.md#workflow-example)).
- Add `NIGHT` as a grouping keyword, Pre-processing only, so each night's
  lights get that night's flats.
- CosmeticCorrection on, PSF Signal Weight, one registration reference for
  all filters, LocalNormalization, ESD rejection, autocrop.

WBPP integrates each exposure time separately, so if a filter was shot at more
than one exposure (here H and S at 300 s and 600 s) you get one master per
exposure. Step 2 merges them.

Save a project (File → Save Project) now and again from here on. It keeps
the images, their history and your process icons together if PixInsight
crashes.

## 2. One master per filter

When WBPP registered every group to the same reference image, the masters
are already aligned; check the log for a single "Reference image" line. If
they aren't, run StarAlignment first.

**Best: integrate all the subs together.** More frames give better rejection
than two masters.

1. In ImageIntegration, add the `*_r.xisf` files from each exposure's folder
   under `registered\`.
2. Weights: PSF Signal Weight. This gives the longer subs more weight without
   any arithmetic.
3. Normalization: Additive with scaling. Rejection: generalized ESD or
   Winsorized sigma clipping.
4. WBPP's LocalNormalization used a different reference for each group, so
   either turn LN off here or run LocalNormalization again on all the frames
   against one reference. Don't mix the groups' `.xnml` files.

**Quick: integrate the two masters.** ImageIntegration with the two master
files, no rejection, PSF Signal Weight, Additive with scaling. In PixelMath,
weight by total exposure rather than 50/50; for 40 min and 120 min that is
`0.25*H300 + 0.75*H600`.

## 3. Shared steps, on all six masters

1. **DynamicCrop.** Set it up on one master and drag its triangle onto the
   others, so all six get the same crop. Cut off the ragged edges where
   nights didn't overlap.
2. **BlurXTerminator**, Correct Only.
3. **ImageSolver** on each master.
4. **SpectrophotometricFluxCalibration (SPFC).** For H, S and O, enter the
   filter's centre wavelength and bandwidth. Wrong values here put that
   channel on the wrong brightness scale.
5. **MultiscaleGradientCorrection (MGC).** See below.

### MGC and Gradient Scale

MGC compares your image with the MARS sky survey, which is free of gradients.
SPFC put your image in the same flux units so the two can be compared.

Gradient Scale is the size, in pixels, that MGC treats as gradient. Anything
larger is taken from MARS and your version of it is replaced; anything smaller
is kept from your image.

On this rig 1024 px is about 22′ and 4096 px about 88′. The Heart Nebula
fills most of the frame, so at 1024 its broad glow counted as gradient and
was flattened or left blotchy. At **4096** only broad trends are removed, such
as a corner-to-corner slope or a moonlit wash, and the nebula's large-scale
shape stays. The cost is that a gradient smaller than about 1.5° stays too.

To check the setting:

- Turn on MGC's option to output the gradient model. It should be a smooth
  plane or gentle curve with no trace of the nebula.
- Look at the result with a hard STF. If one side is still brighter, try
  2048, or follow MGC with a gentle DBE with samples only on true background.

If MGC won't take a master, use GradientCorrection, DBE or GraXpert on it.
Go carefully on H: the faint outer nebula is easy to mistake for gradient.

## 4. Nebula: H, S, O

1. **BlurXTerminator**, normal mode. Sharpen Stars about 0.25, Sharpen
   Nonstellar about 0.4–0.5. These stars are thrown away, so the nonstellar
   setting is the one that matters.
2. **NoiseXTerminator.** Light on H (0.3–0.5), a bit more on S, about 0.5–0.6
   on O. Weaker channels take more.
3. **StarXTerminator** with "Generate star image" off. Keep the starless H,
   S and O.
4. **Stretch each starless master on its own.** Use MultiscaleAdaptiveStretch,
   GHS or HistogramTransformation, one image at a time. Aim for:
   - about the same background in all three, roughly 0.08–0.12 (check with
     the readout or Statistics);
   - each channel's own nebula detail clearly visible.

   O will look fainter and noisier than H. That's normal.
5. **ChannelCombination** (RGB) of the stretched images. For SHO: R = S,
   G = H, B = O. For HOO: R = H, G = O, B = O.
6. **NarrowbandNormalization** (Process → ColorCalibration, a process since
   PixInsight 1.8.9-2). It works on one combined RGB image, preferably
   starless and stretched. Set the palette to match step 5 and turn on
   real-time preview. Adjust the OIII and SII amounts first, then shadows,
   highlights and brightness. Its SCNR setting controls green. Check its
   documentation (the icon at the bottom of the process window) for the
   channel order it expects; if the preview is badly off-colour from the
   start, that is the first thing to check.
7. **Contrast:** LocalHistogramEqualization (below), then small curves and
   saturation changes.

### Why not LinearFit, or one stretch on the combined image

Here O's median was about 10× that of H and S: longer subs (1200 s), moonlight
(brightest near 500 nm, right where OIII is), a WBPP pedestal and the SPFC
scaling all add up. It isn't a fault in the data.

- **Combining linear and stretching once** turned the whole image bright
  blue, because O's higher background dominated.
- **LinearFit with H as reference** fixed the background but also scaled O's
  real signal down, because LinearFit fits every pixel and the nebula is far
  brighter in Ha than in OIII. The blue all but disappeared from the final
  image.

What differs between the channels is mainly the background, an added
offset, and stretching each channel separately (step 4) deals with that and
with O's weaker signal at once. If you want to stay linear until after
combining, subtract O's offset without scaling it: PixelMath on O with
`$T - med($T) + med(H)`, where `H` is the starless H view. A single stretch
still keeps O faint, though.

### LocalHistogramEqualization

The defaults look rippled because Amount 1.0 is full strength, the 8-bit
histogram splits smooth gradients into contour lines, and with nothing masked
the background noise gets boosted too. A zoomed-out real-time preview also
runs on a reduced copy of the image, so it doesn't match the result.

1. Work on the stretched, starless nebula. On a colour image LHE changes only
   lightness.
2. Make a luminance mask: extract the lightness (ChannelExtraction, CIE L*),
   stretch it a bit more, blur it 3–5 px with Convolution, and apply it to the
   image. The nebula is exposed and the dark background protected.
3. Draw a Preview over the core and a lobe edge and look at it at 1:1.
4. Settings:

   | setting | default | use |
   |---|---|---|
   | Kernel Radius | 64 | 150–250 for the large structures |
   | Contrast Limit | 2.0 | 1.5 |
   | Amount | 1.0 | 0.2–0.3 |
   | Histogram Resolution | 8-bit | 10-bit or 12-bit |
   | Circular Kernel | on | on |

5. For more depth, use two light passes rather than one strong one: radius
   about 250 at amount 0.25, then radius 48–64 at 0.1–0.15.

Dark rings around bright areas mean the Contrast Limit or Amount is too high.
A blotchy background means the mask isn't dark enough there; stretch it
harder.

## 5. Stars: R, G, B

1. **ChannelCombination** of R, G and B, still linear.
2. **ImageSolver**, then **SpectrophotometricColorCalibration (SPCC)**. Sensor
   curve: Sony IMX571. Filters: your R, G and B filters.
3. **BlurXTerminator**, normal mode, on the linear image. This is where star
   size is best controlled:
   - Sharpen Stars 0.4–0.6 for tighter cores;
   - Adjust Star Halos −0.2 to −0.5 to shrink the halos of bright stars
     without touching faint ones.
4. **StarXTerminator** with "Generate star image" on. Keep the stars image
   and discard the starless one.
5. **Stretch gently.** Bright stars swell during the stretch.
   - ArcsinhStretch keeps colour in the cores; keep the stretch factor low.
   - With GHS, a high local intensity (b about 6–12) with the symmetry point
     near the faint end lifts faint stars and holds back bright ones.

   Add a little saturation with CurvesTransformation.

### Reducing large stars

MorphologicalTransformation erosion shrinks every star by the same number of
pixels, so small stars vanish while large ones barely change. Large stars
look large because of their halos, so do steps 3 and 5 above first, and dim
the stars when merging (step 6). Often that is enough.

If it isn't, aim MT at the large stars only:

1. StarMask on the stars image with Scale 6–7 and Small-scale 2–3, so small
   stars are left out. Grow the mask a little so it covers the halos.
2. Apply the mask to the stars image.
3. MorphologicalTransformation: Morphological Selection at 0.2–0.35 (gentler
   than Erosion, fewer dark rings), circular structuring element 7×7 or 9×9,
   Amount 0.5–0.7, 1–2 iterations.

Bill Blanshan's Star Reduction script (installed from his PixInsight
repository) has a "Halo" method for the same job.

## 6. Merge and finish

1. Open PixelMath with the nebula as the target and screen the stars on,
   where `stars` is the stars image's view name:
   ```
   ~((~$T)*(~stars))
   ```
   To dim the stars, scale them in the same expression; 0.6–0.8 works:
   ```
   ~((~$T)*(~(stars*0.7)))
   ```
2. A light NoiseXTerminator, then small curves and saturation tweaks.
3. Save the working file as XISF, an archive copy as 16-bit TIFF and a JPG
   to share.

## Process icons

Icons are in `D:\Astro\Pixinsight\Icons.xpsm` (Process → Process Icons → Load
Process Icons). Once you've tuned them, save them to a new file (Save Process
Icons) so the base set stays as it is. After a crash,
[pi-swap-recover](../pi-swap-recover) can recover images from the swap files.
