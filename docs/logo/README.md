<img src="logo.svg" width="96" align="right" alt="astro-tools logo">

# astro-tools logo

A lowercase *a* whose bowl is a telescope aperture with a bright star in it,
circled by an orbit.

| file | use |
|---|---|
| `logo.svg` | the main logo: mark on a dark rounded tile; works on light and dark pages |
| `logo-small.svg` | simpler tile for 48 px and below (no diffraction spikes, fewer stars) |
| `logo-mark-dark.svg` | mark alone, transparent, for dark backgrounds |
| `logo-mark-light.svg` | mark alone, transparent, for light backgrounds |
| `logo-16.png` … `logo-1024.png` | the tile at 16, 32, 48, 64, 128, 256, 512 and 1024 px (16–48 from `logo-small.svg`) |
| `logo-mark-dark-512.png`, `logo-mark-light-512.png` | transparent marks at 512 px |
| `favicon.ico` | 16–256 px |

Colours: night `#0B1020`, starlight `#E8ECF4`, H-alpha `#FF5A4E`
(`#E5483D` on light backgrounds).

To change the logo, edit `make_logo.py` and run it; it rewrites every file here:

```
pip install resvg-py pillow
python docs/logo/make_logo.py
```
