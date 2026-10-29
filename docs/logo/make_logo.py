"""Regenerate the astro-tools logo files in this folder.

Writes the SVGs from the geometry below, then renders the PNGs and favicon.ico from them.
  pip install resvg-py pillow
  python docs/logo/make_logo.py
"""
import io, os
import resvg_py
from PIL import Image

OUT = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)

INK = dict(fg="#E8ECF4", dim="#8C97B3", accent="#FF5A4E", core="#FFF4EE")     # for dark backgrounds
INK_LIGHT = dict(fg="#141B30", dim="#5B6684", accent="#E5483D", core="#FFF4EE")  # for light backgrounds
ORBIT = 'cx="96" cy="112" rx="88" ry="26" transform="rotate(-18 96 112)"'


def star(cx, cy, r):
    p = r * 0.12
    f = lambda v: f"{v:g}"
    return (f"M{f(cx)} {f(cy - r)} Q{f(cx + p)} {f(cy - p)} {f(cx + r)} {f(cy)} Q{f(cx + p)} {f(cy + p)} {f(cx)} {f(cy + r)} "
            f"Q{f(cx - p)} {f(cy + p)} {f(cx - r)} {f(cy)} Q{f(cx - p)} {f(cy - p)} {f(cx)} {f(cy - r)} Z")


def mark(c, small=False, uid="a"):
    """The monogram in its 200x200 box. Orbit gaps are masks, so it works on any background."""
    ow = 6 if small else 4
    sr = 30 if small else 26
    s = f'''<defs>
<clipPath id="{uid}-front"><rect x="-60" y="112" width="320" height="120" transform="rotate(-18 96 112)"/></clipPath>
<mask id="{uid}-back" maskUnits="userSpaceOnUse" x="0" y="0" width="200" height="200">
<rect width="200" height="200" fill="#fff"/>
<circle cx="92" cy="112" r="52" fill="none" stroke="#000" stroke-width="30"/>
<rect x="129" y="50" width="30" height="122" rx="15" fill="#000"/>
</mask>
<mask id="{uid}-cut" maskUnits="userSpaceOnUse" x="0" y="0" width="200" height="200">
<rect width="200" height="200" fill="#fff"/>
<g clip-path="url(#{uid}-front)"><ellipse {ORBIT} fill="none" stroke="#000" stroke-width="{ow + 8}"/></g>
</mask>
</defs>
<ellipse {ORBIT} fill="none" stroke="{c['accent']}" stroke-width="{ow}" opacity="0.55" mask="url(#{uid}-back)"/>
<g mask="url(#{uid}-cut)">
<circle cx="92" cy="112" r="52" fill="none" stroke="{c['fg']}" stroke-width="18"/>
<rect x="135" y="56" width="18" height="110" rx="9" fill="{c['fg']}"/>
</g>
<ellipse {ORBIT} fill="none" stroke="{c['accent']}" stroke-width="{ow}" clip-path="url(#{uid}-front)"/>
'''
    if not small:
        s += f'<path d="M92 72 V152 M52 112 H132" fill="none" stroke="{c["accent"]}" stroke-width="1.5" stroke-linecap="round" opacity="0.7"/>\n'
    s += f'<path d="{star(92, 112, sr)}" fill="{c["accent"]}"/>\n<circle cx="92" cy="112" r="{5 if small else 4}" fill="{c["core"]}"/>\n'
    if not small:
        s += (f'<path d="{star(40, 42, 9)}" fill="{c["fg"]}"/>\n'
              f'<path d="{star(168, 30, 6)}" fill="{c["dim"]}"/>\n'
              f'<path d="{star(30, 172, 5)}" fill="{c["fg"]}" opacity="0.75"/>\n')
    else:
        s += f'<path d="{star(40, 42, 12)}" fill="{c["fg"]}"/>\n'
    return s


def svg(body, label):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 200" role="img" aria-label="{label}">\n'
            f'<title>{label}</title>\n{body}</svg>\n')


def tile(small=False):
    # dark rounded square, mark scaled to sit centred with breathing room; reads on light and dark pages
    inner = mark(INK, small, uid="s" if small else "t")
    k = 0.92 if small else 0.86
    return (f'<rect width="200" height="200" rx="44" fill="#0B1020"/>\n'
            f'<g transform="translate(100 100) scale({k}) translate(-98 -101)">\n{inner}</g>\n')


files = {
    "logo.svg": svg(tile(), "astro-tools"),
    "logo-small.svg": svg(tile(small=True), "astro-tools"),
    "logo-mark-dark.svg": svg(mark(INK, uid="d"), "astro-tools"),
    "logo-mark-light.svg": svg(mark(INK_LIGHT, uid="l"), "astro-tools"),
}
for name, text in files.items():
    open(os.path.join(OUT, name), "w").write(text)


def png(name, size):
    data = bytes(resvg_py.svg_to_bytes(svg_string=files[name], width=size, height=size))
    return Image.open(io.BytesIO(data)).convert("RGBA")


icons = []
for size in (16, 32, 48, 64, 128, 256, 512, 1024):
    im = png("logo-small.svg" if size <= 48 else "logo.svg", size)
    im.save(os.path.join(OUT, f"logo-{size}.png"), optimize=True)
    if size <= 256:
        icons.append(im)
icons[-1].save(os.path.join(OUT, "favicon.ico"), sizes=[(i.width, i.height) for i in icons], append_images=icons[:-1])
for name in ("logo-mark-dark.svg", "logo-mark-light.svg"):
    png(name, 512).save(os.path.join(OUT, name.replace(".svg", "-512.png")), optimize=True)
print("\n".join(f"{f:28s} {os.path.getsize(os.path.join(OUT, f)):>8d}" for f in sorted(os.listdir(OUT))))
