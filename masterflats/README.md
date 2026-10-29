<a href="../README.md"><img src="../docs/logo/logo.svg" width="64" align="right" alt="astro-tools"></a>

# masterflats

Turns each night's raw flats into one master flat per filter, written next to
the raw flats, so you don't have to keep hundreds of flats around. It is safe
to run every day: masters that are already up to date are left alone.

```
D:\Astro\Downloads\Heart Nebula\NIGHT_2026-10-27\FLAT\H\
    2026-10-28_06-28-24_H_-10.00_10.00s_0000.xisf    25 raw flats...
    MasterFlat_2026-10-27_H_10s.xisf                 ...become this
```

## Install

Same as the rest of [astro-tools](../README.md): clone the repo, make a venv and
`pip install -r requirements.txt`. masterflats uses the file readers in
`../subselect`, so keep the two folders side by side.

## Use

```
# see what it would do; changes nothing
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master" -n

# build any new or changed masters
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master"

# and delete the raw flats once they are in a master
python masterflats/masterflats.py "D:\Astro\Downloads" --darks "D:\Astro\Darks -10C\master" --delete-raw
```

Linux, macOS and WSL take the same arguments with `/` paths. Quote paths with
spaces. Wildcards (`"D:\Astro\Downloads\*\NIGHT_2026-10-2*"`) work in cmd and
PowerShell too. Folders whose name starts with `_selected` (subselect's
output) are skipped.

Options:

| option | |
|---|---|
| `--darks DIR` | folder of master darks / biases (IMAGETYP `Master Dark` / `Master Bias`); repeatable |
| `--delete-raw` | delete raw flats once they are in a master |
| `-n`, `--dry-run` | report only |
| `--min-flats N` | fewest flats to build a master from (default 3) |
| `--kappa K` | sigma-clipping threshold (default 3) |
| `-j N` | threads reading flats (default 4) |

## How a master is made

1. **Grouping.** Raw flats (IMAGETYP `FLAT`) are grouped by folder, night,
   filter, exposure, binning, gain, offset and rotator angle. Night comes from a
   `NIGHT_YYYY-MM-DD` folder in the path, or local time minus 12 h. Groups with
   fewer than `--min-flats` flats are skipped.
2. **Flat-dark.** The camera offset has to come off before flats are divided
   into lights, so each group gets one of these, in order:
   - flat-darks of the same night, exposure and gain anywhere in the tree
     (IMAGETYP `DARK` or `DARKFLAT`; several are median-combined);
   - a master dark of the same exposure from `--darks`;
   - a master bias from `--darks`;
   - the master dark nearest in exposure. On a cooled CMOS the dark current
     in a 10 s flat is negligible.

   Masters must be within 2 °C of the flats' set temperature (from `SET-TEMP`
   or a folder name like `-10C`). With nothing suitable, the group is **not
   built**.
3. **Combine.** Each flat minus the dark is divided by its median. Every pixel
   is then a mean of the frames within `kappa` × MAD-sigma of their median.
   That removes cosmic rays, satellites and hot pixels. The result is scaled to
   the mean flat level.
4. **Output.** Float32 XISF in [0, 1], zstd-compressed with byte shuffling
   (the same as PixInsight's "zstd + byte shuffling" save option; about 60 MB
   instead of 104 MB for a 26 MP frame, lossless), with IMAGETYP
   `Master Flat`, imageType `MasterFlat`, `FILTER`, `EXPTIME`, `NCOMBINE`, and
   camera, rotator and date keywords copied from the first flat. WBPP uses it
   as-is. It is written to a temp file and then renamed, so an interrupted run
   never leaves half a master.

Checked against WBPP 3.1.2 on the same 25 flats and 10 s master dark (ESD
rejection, multiplicative normalisation): the large-scale shape agrees to
0.001%, single pixels to 0.06% RMS (worst pixel 1%), from the slightly
different outlier rejection.

## Idempotency

Each master records the names and sizes of the flats and flat-darks that went
into it. On the next run, a group is:

- **up to date** when the inputs match. That includes when some or all raw flats
  were deleted after the build: the master is never rebuilt from fewer flats;
- **compressed in place** (not rebuilt) when it is up to date but was written
  uncompressed by an older version;
- **rebuilt** when it has flats the master doesn't (more flats taken that
  morning), or when a different flat-dark now applies (e.g. you added flat-darks).

Masters made by other software (WBPP etc.) are never touched.

## With subselect

`subselect.py --copy` copies a night's master flat into
`_selected/NIGHT_<date>/FLAT/` instead of its raw flats. It also tells you if
raw flats from an earlier run are still sitting there.

## Tests

```
pip install pytest
python -m pytest masterflats
```
