# pi-swap-recover

Gets your images back from PixInsight's swap files after a crash.

PixInsight saves each image's undo history as `~PI~<id>~<view>.swp` files, plus
`-partNNN.swp` files, in its swap folders. They're normally deleted on exit,
but they survive a crash. This script converts every history state it finds
into a normal XISF file.

## Where the swap files are

PixInsight → **Edit → Global Preferences → Directories and Network** → swap
directories. The default is the system temp folder: `%TEMP%` on Windows
(usually `C:\Users\<you>\AppData\Local\Temp`), `/tmp` on Linux and macOS.
The setting is stored in `%APPDATA%\Pleiades\PixInsight.ini` as
`SwapDirectory001=D:/Temp` etc.; if you don't give a swap folder, the script
reads it from there.

If you set more than one swap directory, the `-partNNN` files are spread over
all of them. Copy (or pass) all of them.

## Copy the swap files first

PixInsight may delete or reuse them once it restarts, and a running PixInsight
keeps them locked. So before anything else, with PixInsight closed:

```
mkdir D:\PI_recovery\swap
copy "D:\Temp\~PI~*.swp" D:\PI_recovery\swap\
```

Use your own swap folder in place of `D:\Temp` (`copy "%TEMP%\~PI~*.swp" ...`
for the default). `copy` keeps the file times, which the script uses to order
the states. Or let the script do it with `--copy-to D:\PI_recovery\swap`.

## Install (Windows)

Needs Python 3 (from python.org or the Microsoft Store). In a Command Prompt,
in the folder with `pi_swap_recover.py`:

```
py -m venv .venv
.venv\Scripts\activate
pip install numpy pillow
```

On Linux/macOS: `python3 -m venv .venv && . .venv/bin/activate && pip install numpy pillow`.

## Run

```
python pi_swap_recover.py D:\PI_recovery\swap D:\PI_recovery\recovered --previews
python pi_swap_recover.py SWAP OUT --latest              # only the newest state of each image
python pi_swap_recover.py SWAP OUT --view MasterLight_H  # one image's history
python pi_swap_recover.py D:\Swap1 E:\Swap2 OUT          # several swap folders
python pi_swap_recover.py OUT --copy-to D:\PI_recovery\swap --previews
                                       # swap folders from PixInsight's settings, copied first
```

`D:\Temp`, `D:/Temp` and `"C:\My Swap"` (quotes for spaces) all work.

Output files are named `<date>_<time>_<view>_<id>.xisf`. `--previews` adds an
auto-stretched PNG next to each one, so you can find the one you want without
opening them all. Requires numpy, plus pillow for `--previews`.

States that can't be read (a part missing, a file still locked by PixInsight,
not a swap file) are listed as SKIPPED and the rest are still recovered.
Swap compression (Global Preferences) isn't supported; such states show up as
truncated.

## What you get back

- Each state is the image **before** a process was applied to it. So the
  newest state of a view is one step behind what was on screen when PixInsight
  crashed. Work that was never followed by another process (the image on
  screen) isn't in the swap files.
- **Pixels only.** The FITS keywords, astrometric solution and processing
  history aren't in the swap files. Re-run ImageSolver if you need the solution
  back.

## Tests

```
pip install pytest
python -m pytest pi-swap-recover -q
```

## File format (as reverse-engineered from PixInsight 1.9.5)

- **Header, 128 bytes in the main file, little-endian:** magic `egamiwar`.
  - byte 12: float flag
  - byte 15: bits per sample
  - offset 0x1C: width, height, channels, colour space (1 = RGB), as uint32
  - offset 0x2C: RGB working space values
  - offset 0x68: number of parts (uint32)
- **Pixel data:** split evenly across the main file and `part001`…`partNNN`,
  one per swap directory. Each part holds its slice of **every channel in
  turn**: part k = [channel 0 slice k, channel 1 slice k, …].
