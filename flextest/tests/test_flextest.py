"""Tests for flextest.py / skybg.py on synthetic XISF subs. Run: python -m pytest flextest -q"""
import os, shutil, struct, subprocess, sys
from datetime import datetime, timedelta
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
import flextest as ft

N = 512                  # image size, px
DRIFT = (6.0, 3.0)       # px per sub, also the trail inside each sub (no gaps between subs)
EXP = 60.0
T0 = datetime(2026, 10, 28, 1, 0, 0)


# ---------------------------------------------------------------- synthetic data
def write_xisf(path, img, kw, codec=None, fmt="UInt16"):
    if fmt == "UInt16":
        raw = np.clip(np.round(img), 0, 65535).astype("<u2").tobytes(); isz = 2
    else:
        raw = (img / 65535).astype("<f4").tobytes(); isz = 4
    data, comp = raw, ""
    if codec:
        base, sh = codec.split("+")[0], codec.endswith("+sh")
        if sh:
            raw2 = np.frombuffer(raw, np.uint8).reshape(-1, isz).T.tobytes()
        else:
            raw2 = raw
        if base.startswith("lz4"):
            import lz4.block
            data = lz4.block.compress(raw2, store_size=False, mode="high_compression" if base == "lz4hc" else "default")
        elif base == "zstd":
            import zstandard
            data = zstandard.ZstdCompressor().compress(raw2)
        else:
            import zlib
            data = zlib.compress(raw2)
        comp = f' compression="{codec}:{len(raw)}' + (f':{isz}"' if sh else '"')
    keys = "".join(f'<FITSKeyword name="{k}" value="{v}" comment=""/>' for k, v in kw.items())
    off = 4096
    xml = (f'<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0"><Image geometry="{img.shape[1]}:{img.shape[0]}:1" '
           f'sampleFormat="{fmt}" colorSpace="Gray" location="attachment:{off}:{len(data)}"{comp}>{keys}</Image></xisf>').encode()
    assert len(xml) + 16 <= off
    with open(path, "wb") as f:
        f.write(b"XISF0100" + struct.pack("<II", len(xml), 0) + xml)
        f.write(b"\0" * (off - 16 - len(xml)) + data)


def star_field(n=90, seed=1):
    rng = np.random.default_rng(seed)
    pts = []
    while len(pts) < n:
        p = rng.uniform(45, N - 45 - 30, 2)        # room for the accumulated drift
        if all(np.hypot(*(p - q)) > 30 for q in pts):
            pts.append(p)
    return np.array(pts), rng.uniform(400, 4000, n)


def render(pts, amp, offset, trail, sigma=1.4, seed=0):
    """Gaussian stars, each smeared along `trail` starting at `offset`."""
    rng = np.random.default_rng(seed)
    img = 1000 + rng.normal(0, 10, (N, N))
    yy, xx = np.mgrid[-20:21, -20:21].astype(float)
    steps = np.linspace(0, 1, 25)
    for (x, y), a in zip(pts + offset, amp):
        ix, iy = int(round(x)), int(round(y))
        fx, fy = x - ix, y - iy
        cut = sum(np.exp(-((xx - fx - s * trail[0]) ** 2 + (yy - fy - s * trail[1]) ** 2) / (2 * sigma ** 2)) for s in steps)
        img[iy - 20:iy + 21, ix - 20:ix + 21] += a * cut / len(steps)
    return img


def header(i, filt="H"):
    t = T0 + timedelta(seconds=i * EXP)
    return {"FILTER": f"'{filt}'", "EXPTIME": EXP, "EGAIN": 0.25, "XPIXSZ": 3.76, "FOCALLEN": 600.0,
            "DATE-OBS": f"'{t.isoformat(timespec='microseconds')}0'", "DATE-LOC": f"'{(t - timedelta(hours=5)).isoformat()}'",
            "PIERSIDE": "'West'", "FOCPOS": 1000, "FOCTEMP": 5.0}


@pytest.fixture(scope="session")
def subs(tmp_path_factory):
    """4 trailed+drifting subs and 4 round, static subs, in dirs with spaces in the name."""
    root = tmp_path_factory.mktemp("data")
    pts, amp = star_field()
    out = {"trail": [], "round": []}
    codecs = [None, "lz4+sh", "zstd+sh", "lz4hc+sh"]
    for kind, d in [("trail", np.array(DRIFT)), ("round", np.zeros(2))]:
        dd = root / f"{kind} subs"; dd.mkdir()
        for i in range(4):
            j = i if kind == "trail" else i + 10
            img = render(pts, amp, i * d, d, seed=j)
            p = dd / f"2026-10-28_{kind}_{i:04d}.xisf"
            write_xisf(p, img, header(j), codec=codecs[i])
            out[kind].append(str(p))
    dark = root / "dark.xisf"
    write_xisf(dark, np.full((N, N), 1000.0) + np.random.default_rng(5).normal(0, 3, (N, N)), header(0, "D"))
    out["dark"] = str(dark)
    out["root"] = str(root)
    return out


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    monkeypatch.setattr(ft, "CACHE", str(tmp_path / "cache"))
    monkeypatch.setenv("FLEXTEST_CACHE", str(tmp_path / "cache"))
    return tmp_path / "cache"


# ---------------------------------------------------------------- unit tests
@pytest.mark.parametrize("codec", [None, "lz4", "lz4+sh", "lz4hc+sh", "zstd", "zstd+sh", "zlib+sh"])
@pytest.mark.parametrize("fmt", ["UInt16", "Float32"])
def test_read_xisf_codecs(tmp_path, codec, fmt):
    img = np.random.default_rng(0).integers(0, 60000, (37, 53)).astype(np.float32)
    p = tmp_path / "x.xisf"
    write_xisf(p, img, {"FILTER": "'Ha &amp; O3'", "EXPTIME": 300.0}, codec=codec, fmt=fmt)
    got, hdr = ft.read_xisf(str(p))
    assert got.shape == (37, 53)
    np.testing.assert_allclose(got, img, atol=0.05 if fmt == "Float32" else 0)
    assert hdr["FILTER"] == "Ha & O3" and float(hdr["EXPTIME"]) == 300


def test_read_xisf_rejects_non_xisf(tmp_path):
    p = tmp_path / "x.fits"; p.write_bytes(b"SIMPLE  =" + b" " * 100)
    with pytest.raises(ValueError):
        ft.read_xisf(str(p))


def test_adaptive_moments_round_and_trailed():
    yy, xx = np.mgrid[:25, :25].astype(float)
    g = np.exp(-((xx - 12.3) ** 2 + (yy - 11.8) ** 2) / (2 * 1.5 ** 2))
    cx, cy, C = ft.adaptive_moments(g, 12, 12)
    assert abs(cx - 12.3) < 0.01 and abs(cy - 11.8) < 0.01
    np.testing.assert_allclose(C, np.diag([2.25, 2.25]), atol=0.01)
    # line of length L along x adds L^2/12 to Cxx
    L = 6.0
    t = sum(np.exp(-((xx - 12 - s) ** 2 + (yy - 12) ** 2) / (2 * 1.5 ** 2)) for s in np.linspace(-L / 2, L / 2, 61))
    _, _, C = ft.adaptive_moments(t, 12, 12)
    assert 0.9 < (C[0, 0] - C[1, 1]) / (L ** 2 / 12) < 1.5      # Gaussian weight over-reads a flat-topped trail a bit
    assert abs(C[0, 1]) < 0.01


def test_measure_trail(subs):
    s, hdr, shape = ft.measure(subs["trail"][1])
    assert shape == (N, N) and len(s) > 60 and hdr["FILTER"] == "H"
    st = ft.shape_stats(s)
    L = np.hypot(*DRIFT)
    assert st["trail"] == pytest.approx(L, rel=0.3)
    assert st["pa"] == pytest.approx(np.degrees(np.arctan2(DRIFT[1], DRIFT[0])), abs=5)
    r = ft.shape_stats(ft.measure(subs["round"][1])[0])
    assert r["trail"] < 1.5 and r["ecc"] < st["ecc"]
    assert r["fwhm"] == pytest.approx(2.3548 * 1.4, rel=0.15)


def test_drift_recovered(subs):
    m = [ft.measure(f) for f in subs["trail"]]
    for i in range(1, len(m)):
        r = ft.register(m[i - 1][0], m[i][0], m[i][2])
        assert r is not None and r["nmatch"] > 40
        assert r["dx"] == pytest.approx(DRIFT[0], abs=0.1) and r["dy"] == pytest.approx(DRIFT[1], abs=0.1)
        assert abs(r["rot"]) < 0.5
    r = ft.register(m[0][0], m[3][0], m[3][2])
    assert (r["dx"], r["dy"]) == pytest.approx((3 * DRIFT[0], 3 * DRIFT[1]), abs=0.15)
    a, b = (ft.measure(f)[0] for f in subs["round"][:2])
    r = ft.register(a, b, (N, N))
    assert abs(r["dx"]) < 0.05 and abs(r["dy"]) < 0.05


def test_cache_hit_and_key_ignores_path(subs, cache, tmp_path):
    f = subs["trail"][0]
    s1, h1, sh1 = ft.measure(f)
    files = os.listdir(cache)
    assert len(files) == 1 and files[0].endswith(".npz")
    cp = tmp_path / "elsewhere" / "copy.xisf"; cp.parent.mkdir()
    shutil.copy(f, cp)
    assert ft.cache_key(str(cp)) == ft.cache_key(f)
    s2, h2, sh2 = ft.measure(str(cp))
    assert len(os.listdir(cache)) == 1
    np.testing.assert_array_equal(s1, s2); assert h1 == h2 and sh1 == sh2
    (cache / files[0]).write_bytes(b"garbage")          # corrupt entry: re-measured, not fatal
    np.testing.assert_array_equal(ft.measure(f)[0], s1)


def test_unwritable_cache_warns(subs, tmp_path, monkeypatch, capsys):
    blocker = tmp_path / "notadir"; blocker.write_text("x")
    monkeypatch.setattr(ft, "CACHE", str(blocker / "cache"))
    monkeypatch.setattr(ft, "_cache_warned", False)
    s, _, _ = ft.measure(subs["round"][0])
    assert len(s) > 60
    assert "can't write cache" in capsys.readouterr().err


def test_expand(tmp_path):
    for n in ["a.xisf", "b.xisf", "c[1].xisf", "d.fits"]:
        (tmp_path / n).write_bytes(b"")
    got = ft.expand([str(tmp_path / "*.xisf")])
    assert [os.path.basename(g) for g in got] == ["a.xisf", "b.xisf", "c[1].xisf"]
    assert ft.expand([str(tmp_path / "c[1].xisf")]) == [str(tmp_path / "c[1].xisf")]   # literal name wins
    assert ft.expand(["no_such_plain_file.xisf"]) == ["no_such_plain_file.xisf"]
    with pytest.raises(SystemExit):
        ft.expand([str(tmp_path / "*.nope")])


# ---------------------------------------------------------------- end-to-end
def run(script, *args, env=None):
    e = dict(os.environ, FLEXTEST_CACHE=env.pop("cache"), **env) if env else dict(os.environ)
    return subprocess.run([sys.executable, os.path.join(SRC, script), *args], capture_output=True, env=e, timeout=300)


def test_e2e_flextest_glob(subs, cache):
    pat = os.path.join(subs["root"], "trail subs", "*.xisf")
    p = run("flextest.py", "--map", pat, env={"cache": str(cache), "PYTHONIOENCODING": "utf-8"})
    out = p.stdout.decode("utf-8")
    assert p.returncode == 0, p.stderr.decode("utf-8", "replace")
    assert out.count("\n 0->1 ") + out.count("\n 1->2 ") + out.count("\n 2->3 ") == 3
    rate = np.hypot(*DRIFT) * ft.scale_of(header(0))        # px/sub, subs are 1 min apart
    line = next(l for l in out.splitlines() if l.strip().startswith("-: 4 subs"))
    assert float(line.split("drift ")[1].split()[0]) == pytest.approx(rate, rel=0.03)
    assert len(os.listdir(cache)) == 4


def test_e2e_off_on_globs(subs, cache):
    p = run("flextest.py", "--off", os.path.join(subs["root"], "round subs", "*.xisf"),
            "--on", os.path.join(subs["root"], "trail*", "*_000[0-3].xisf"),
            env={"cache": str(cache), "PYTHONIOENCODING": "utf-8"})
    assert p.returncode == 0, p.stderr.decode("utf-8", "replace")
    out = p.stdout.decode()
    assert "compensator:" in out and "off: 4 subs" in out and "on: 4 subs" in out
    assert p.stderr.decode().count("measuring") == 8


def test_e2e_no_match_is_clean_error(subs, cache):
    p = run("flextest.py", os.path.join(subs["root"], "*.nothing"), env={"cache": str(cache)})
    assert p.returncode != 0 and b"no files match" in p.stderr and b"Traceback" not in p.stderr


def test_e2e_skybg_glob(subs):
    p = run("skybg.py", subs["dark"], os.path.join(subs["root"], "round subs", "*.xisf"))
    assert p.returncode == 0, p.stderr.decode("utf-8", "replace")
    lines = p.stdout.decode().splitlines()
    assert len(lines) == 4 and all(l.startswith("2026-10-28_round_") and " H 60s " in l for l in lines)
    sky = float(lines[0].split("sky")[1].split()[0])
    assert abs(sky) < 2                                  # light background == dark level


def test_cp1252_redirect_does_not_crash(subs, cache, tmp_path):
    d = tmp_path / "Herz ″Ω"; d.mkdir()
    for f in subs["trail"][:3]:
        shutil.copy(f, d / ("☆″_" + os.path.basename(f)))
    env = {"PYTHONIOENCODING": "cp1252", "PYTHONUTF8": "0"}
    p = run("flextest.py", str(d / "*.xisf"), env={"cache": str(cache), **env})
    assert p.returncode == 0, p.stderr.decode("cp1252", "replace")
    assert b"Drift fit per group" in p.stdout and b"Traceback" not in p.stderr
    p = run("skybg.py", subs["dark"], str(d / "*.xisf"), env={"cache": str(cache), **env})
    assert p.returncode == 0, p.stderr.decode("cp1252", "replace")
    assert len(p.stdout.splitlines()) == 3 and p.stdout.startswith(b"??_")
