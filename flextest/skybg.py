"""Sky background above dark, in e-, and read-noise swamping. Usage: skybg.py DARK LIGHT..."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flextest import read_xisf
dark, _ = read_xisf(sys.argv[1])
for f in sys.argv[2:]:
    img, h = read_xisf(f)
    eg = float(h.get("EGAIN", 0.2429)); d = img - dark; ts = 256; bg = []
    for y in range(0, d.shape[0] - ts + 1, ts):
        for x in range(0, d.shape[1] - ts + 1, ts):
            t = d[y:y + ts, x:x + ts].ravel()
            for _ in range(5):
                md = np.median(t); sd = 1.4826 * np.median(np.abs(t - md)); t = t[np.abs(t - md) < 3 * sd]
            bg.append(np.median(t))
    bg = np.array(bg); sky = np.median(bg[bg <= np.percentile(bg, 10)]) * eg
    print(f"{f.split('/')[-1][:40]}  {h.get('FILTER')} {float(h['EXPTIME']):.0f}s  sky {sky:5.1f} e-"
          f"  ({sky / eg:5.1f} ADU)  sky/RN^2 {sky / 1.44:4.1f} (RN1.2) {sky / 2.25:4.1f} (RN1.5)"
          f"  eff {sky / (sky + 1.44) * 100:3.0f}-{sky / (sky + 2.25) * 100:3.0f}%  rate {sky / float(h['EXPTIME']) * 60:.2f} e-/min")
