import os, re, subprocess, sys, zlib
from types import SimpleNamespace

import lz4.block
import numpy as np
import pytest
import zstandard

import subselect as ss

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "subselect.py")
KW = {"IMAGETYP": "'LIGHT'", "FILTER": "'H'", "EXPTIME": "300.0", "XPIXSZ": "3.76", "FOCALLEN": "800.0",
      "DATE-LOC": "'2026-10-22T23:10:00.000'", "SET-TEMP": "-10.0", "EGAIN": "0.25"}
NIGHT = "NIGHT_2026-10-22"


# ---------------------------------------------------------------- synthetic data
def star_image(seed=0, sx=1.5, sy=1.5, size=384, bg=1000.0, peak=4000.0):
    """Grid of Gaussian stars on a noisy, slightly sloped sky. sx != sy gives elongated stars."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:size, :size].astype(np.float64)
    img = bg + 40 * xx / size
    for gy in range(24, size - 30, 32):
        for gx in range(24, size - 30, 32):
            cx, cy = gx + rng.uniform(-1, 1), gy + rng.uniform(-1, 1)
            img += peak * np.exp(-0.5 * (((xx - cx) / sx) ** 2 + ((yy - cy) / sy) ** 2))
    img += rng.normal(0, 10, img.shape)
    return np.clip(img, 0, 65535).astype(np.uint16)


def write_xisf(path, img, kw=None, comp=None):
    raw = img.astype("<u2").tobytes()
    data, attr = raw, ""
    if comp:
        src = np.frombuffer(raw, np.uint8).reshape(-1, 2).T.tobytes() if comp.endswith("+sh") else raw
        data = (lz4.block.compress(src, store_size=False) if comp.startswith("lz4") else
                zstandard.ZstdCompressor().compress(src) if comp.startswith("zstd") else zlib.compress(src))
        attr = f' compression="{comp}:{len(raw)}' + (":2" if comp.endswith("+sh") else "") + '"'
    keys = "".join(f'<FITSKeyword name="{k}" value="{v}" comment=""/>' for k, v in (kw or KW).items())
    h, w = img.shape

    def xml(off):
        return (f'<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0"><Image geometry="{w}:{h}:1" '
                f'sampleFormat="UInt16" colorSpace="Gray" location="attachment:{off:010d}:{len(data)}"{attr}>'
                f'{keys}</Image></xisf>').encode("utf-8")
    off = 16 + len(xml(0))
    with open(path, "wb") as f:
        f.write(b"XISF0100" + len(xml(off)).to_bytes(4, "little") + b"\0" * 4 + xml(off) + data)
    return path


def write_fits(path, img, kw=None):
    cards = [("SIMPLE", "T"), ("BITPIX", "16"), ("NAXIS", "2"), ("NAXIS1", str(img.shape[1])),
             ("NAXIS2", str(img.shape[0])), ("BZERO", "32768"), ("BSCALE", "1")] + list((kw or KW).items())
    hdr = "".join(f"{k:<8}= {v:>20} / test".ljust(80) for k, v in cards) + "END".ljust(80)
    hdr = hdr.ljust(-(-len(hdr) // 2880) * 2880)
    data = (img.astype(np.int32) - 32768).astype(">i2").tobytes()
    with open(path, "wb") as f:
        f.write(hdr.encode("ascii") + data + b"\0" * (-len(data) % 2880))
    return path


def kw(**over):
    d = dict(KW)
    d.update({k.replace("_", "-"): v for k, v in over.items()})
    return d


@pytest.fixture
def tree(tmp_path):
    """src/NIGHT_2026-10-22/{LIGHT/H, FLAT}: 4 good lights (raw, lz4+sh, zstd+sh, FITS), 1 trailed, flats."""
    light, flat = tmp_path / "src" / NIGHT / "LIGHT" / "H", tmp_path / "src" / NIGHT / "FLAT"
    light.mkdir(parents=True)
    flat.mkdir(parents=True)
    write_xisf(light / "g1.xisf", star_image(1))
    write_xisf(light / "g2.xisf", star_image(2), comp="lz4+sh")
    write_xisf(light / "g3.xisf", star_image(3), comp="zstd+sh")
    write_fits(light / "g4.fits", star_image(4))
    write_xisf(light / "bad.xisf", star_image(5, sx=3.5, sy=1.2))
    small = np.full((16, 16), 20000, np.uint16)
    fk = kw(IMAGETYP="'FLAT'", EXPTIME="10.0", DATE_LOC="'2026-10-23T06:35:48.000'")
    write_xisf(flat / "flat_H_1.xisf", small, fk)
    write_xisf(flat / "flat_H_2.xisf", small, fk)
    write_xisf(flat / "flat_O_1.xisf", small, dict(fk, FILTER="'O'"))
    write_xisf(flat / "darkflat_1.xisf", small, dict(fk, IMAGETYP="'DARKFLAT'"))
    return tmp_path / "src"


# ---------------------------------------------------------------- readers
@pytest.mark.parametrize("comp", [None, "lz4", "lz4+sh", "zstd", "zstd+sh", "zlib+sh"])
def test_read_xisf(tmp_path, comp):
    img = star_image(0, size=64)
    p = write_xisf(str(tmp_path / "a.xisf"), img, comp=comp)
    out, hdr = ss.read_image(p)
    assert out.dtype == np.float32 and np.array_equal(out, img.astype(np.float32))
    assert hdr == ss.read_header(p)
    assert hdr["FILTER"] == "H" and hdr["IMAGETYP"] == "LIGHT" and float(hdr["EXPTIME"]) == 300
    assert hdr["DATE-LOC"] == "2026-10-22T23:10:00.000"


def test_read_fits(tmp_path):
    img = star_image(0, size=64)
    p = write_fits(str(tmp_path / "a.fits"), img)
    out, hdr = ss.read_image(p)
    assert np.array_equal(out, img.astype(np.float32))
    assert ss.read_header(p)["FILTER"] == "H" and float(hdr["SET-TEMP"]) == -10


def test_read_not_xisf(tmp_path):
    p = tmp_path / "x.xisf"
    p.write_bytes(b"garbage" * 10)
    with pytest.raises(ValueError):
        ss.read_header(str(p))


# ---------------------------------------------------------------- measurement
def test_measure_round_vs_elongated(tmp_path):
    good = ss.measure(write_xisf(str(tmp_path / "good.xisf"), star_image(1)))
    bad = ss.measure(write_xisf(str(tmp_path / "bad.xisf"), star_image(2, sx=3.5, sy=1.2)))
    assert good["n"] >= 100 and bad["n"] >= 100
    assert good["ecc"] < 0.5 and bad["ecc"] > 0.85
    assert 3.0 < good["fwhm"] < 4.2                 # 2.355 * 1.5 px
    assert bad["trail"] > 3 * good["trail"]
    assert abs(bad["pa"] % 180) < 5 or abs(bad["pa"] % 180) > 175   # smeared along x
    assert 900 < good["bg"] < 1100


def test_measure_cache(tmp_path, cache, monkeypatch):
    p = write_xisf(str(tmp_path / "a.xisf"), star_image(1))
    m1 = ss.measure(p)
    monkeypatch.chdir(tmp_path)
    m2 = ss.measure("a.xisf")                       # relative path, same file: same cache entry
    assert m1 == m2 and m1["hdr"]["FILTER"] == "H"
    assert [f for f in os.listdir(cache) if f.endswith(".json")] == os.listdir(cache)
    assert len(os.listdir(cache)) == 1


# ---------------------------------------------------------------- metadata
@pytest.mark.parametrize("path", [r"D:\Astro\Heart Nebula\NIGHT_2026-10-22\LIGHT\H\a.xisf",
                                  "D:/Astro/Heart Nebula/NIGHT_2026-10-22/LIGHT/H/a.xisf",
                                  "/mnt/d/Astro/NIGHT_2026-10-22/LIGHT/a.xisf",
                                  r"\\nas\astro\NIGHT_2026-10-22\FLAT\f.xisf"])
def test_night_from_path(path):
    assert ss.night_of(path, {"DATE-LOC": "2027-01-01T03:00:00"}) == "2026-10-22"


def test_night_from_header():
    assert ss.night_of(r"C:\x\a.xisf", {"DATE-LOC": "2026-10-23T03:00:00.123"}) == "2026-10-22"
    assert ss.night_of("/x/a.xisf", {"DATE-LOC": "2026-10-22T21:00:00"}) == "2026-10-22"
    assert ss.night_of("/x/a.xisf", {}) == "unknown"


def test_pixel_scale_and_temp():
    assert ss.pixel_scale({"XPIXSZ": "3.76", "FOCALLEN": "800"}, None) == pytest.approx(0.9694, abs=1e-3)
    assert ss.pixel_scale({}, None) is None and ss.pixel_scale({}, 1.5) == 1.5
    assert ss.temp_of({"SET-TEMP": "-10"}) == -10
    assert ss.temp_of({}, r"D:\Darks - -10C, Gain 100\m.xisf") == -10
    assert ss.temp_of({}, "/darks/Darks - -15C, Gain 100/m.xisf") == -15


# ---------------------------------------------------------------- grading
LIMITS = SimpleNamespace(min_stars=30, fwhm_good=4.5, fwhm_max=6.0, ecc_good=0.72, ecc_max=0.82,
                         stars_good=0.5, stars_min=0.25, bg_good=1.5, bg_max=2.0, grad_good=2.0, grad_max=3.0)


def sub(**over):
    s = dict(n=200, fwhm=3.5, fwhm_as=3.5, ecc=0.4, star_ratio=1.0, bg_ratio=1.0, grad_ratio=1.0)
    s.update(over)
    return s


@pytest.mark.parametrize("over, tier", [
    ({}, "good"),
    ({"fwhm_as": 5.0}, "borderline"), ({"fwhm_as": 6.5}, "reject"),
    ({"ecc": 0.75}, "borderline"), ({"ecc": 0.9}, "reject"),
    ({"star_ratio": 0.4}, "borderline"), ({"star_ratio": 0.1}, "reject"),
    ({"bg_ratio": 1.7}, "borderline"), ({"bg_ratio": 2.5}, "reject"),
    ({"grad_ratio": 2.5}, "borderline"), ({"grad_ratio": 3.5}, "reject"),
    ({"n": 10}, "reject"), ({"fwhm": None, "fwhm_as": None}, "reject"),
    ({"bg_ratio": None, "grad_ratio": None}, "good"),
    ({"fwhm_as": 5.0, "ecc": 0.9}, "reject"),
])
def test_grade(over, tier):
    t, why = ss.grade(sub(**over), LIMITS)
    assert t == tier and (bool(why) == (tier != "good"))


# ---------------------------------------------------------------- paths
def test_expand_globs(tree):
    light = tree / NIGHT / "LIGHT" / "H"
    assert ss.expand([str(light / "*.xisf")]) == sorted(str(p) for p in light.glob("*.xisf"))
    assert ss.expand([str(tree / "NIGHT_*" / "LIGHT")]) == [str(tree / NIGHT / "LIGHT")]
    fwd = str(tree).replace("\\", "/") + "/NIGHT_*/FLAT"            # forward slashes work on Windows too
    assert [os.path.normpath(p) for p in ss.expand([fwd])] == [str(tree / NIGHT / "FLAT")]
    assert ss.expand([str(light)]) == [str(light)]                  # plain paths pass through
    assert ss.expand([str(tree / "nothing*")]) == []


def test_fix_argv(monkeypatch):
    monkeypatch.setattr(ss.os, "name", "nt")
    assert ss.fix_argv(['D:\\Astro\\Heart Nebula" --copy -j', "2"]) == \
        ["D:\\Astro\\Heart Nebula", "--copy", "-j", "2"]
    assert ss.fix_argv(['D:\\Astro\\Heart Nebula"']) == ["D:\\Astro\\Heart Nebula"]
    assert ss.fix_argv(["D:\\a", "-j", "2"]) == ["D:\\a", "-j", "2"]


def test_find_files_excludes_output(tree):
    out = tree / "_selected"
    (out / NIGHT / "LIGHT").mkdir(parents=True)
    write_xisf(out / NIGHT / "LIGHT" / "copy.xisf", star_image(0, size=32))
    (tree / "notes.txt").write_text("x")
    excl = str(out).upper() if os.name == "nt" else str(out)       # Windows paths compare case-insensitively
    found = ss.find_files([str(tree)], {excl + os.sep})
    assert len(found) == 9 and not any("_selected" in p for p in found)
    assert all(os.path.isabs(p) for p in found)
    assert len(ss.find_files([str(tree)], set())) == 10


# ---------------------------------------------------------------- end to end
def check_layout(out, src):
    light, flat = out / NIGHT / "LIGHT", out / NIGHT / "FLAT"
    assert sorted(os.listdir(light)) == ["g1.xisf", "g2.xisf", "g3.xisf", "g4.fits"]
    assert sorted(os.listdir(flat)) == ["flat_H_1.xisf", "flat_H_2.xisf"]                # O flat unused
    assert os.listdir(out / NIGHT / "DARKFLAT") == ["darkflat_1.xisf"]
    s, d = os.stat(src / NIGHT / "LIGHT" / "H" / "g1.xisf"), os.stat(light / "g1.xisf")
    assert abs(s.st_mtime - d.st_mtime) < 0.01 and s.st_size == d.st_size               # copy2 keeps mtime


def test_main_copy(tree, tmp_path, capsys):
    out = tmp_path / "out" / "_selected"
    ss.main([str(tree), "--copy", "--out", str(out) + os.sep, "-j", "2", "--csv", str(tmp_path / "r.csv")])
    o = capsys.readouterr()
    assert "measuring 5 light subs" in o.err
    tiers = [re.search(r"  (good|BORDERLINE|REJECT)\b", l).group(1) for l in o.out.splitlines()
             if l.startswith("2026-10-22")]
    assert sorted(tiers) == ["REJECT"] + ["good"] * 4
    check_layout(out, tree)
    assert not (tmp_path / "out" / "_selected_borderline").exists()
    assert "copied 7 new file(s) (7 planned)" in o.out
    csv = (tmp_path / "r.csv").read_text(encoding="utf-8-sig").splitlines()
    assert csv[0].startswith("night,time,filter") and len(csv) == 6

    ss.main([str(tree), "--copy", "--out", str(out), "-j", "2"])           # re-run: cached, nothing new
    assert "copied 0 new file(s) (7 planned)" in capsys.readouterr().out


def test_main_default_out_rerun(tree, capsys):
    pattern = str(tree / "NIGHT_*")                   # unexpanded wildcard, as cmd.exe/PowerShell pass it
    ss.main([pattern, "--copy", "-j", "2"])
    capsys.readouterr()
    check_layout(tree / NIGHT / "_selected", tree)
    ss.main([pattern, "--copy", "-j", "2"])           # _selected inside the source tree is not rescanned
    o = capsys.readouterr()
    assert "measuring 5 light subs" in o.err and "copied 0 new file(s)" in o.out


def test_main_wildcard_files(tree, capsys):
    ss.main([str(tree / NIGHT / "LIGHT" / "H" / "g*.xisf"), "-j", "2"])
    o = capsys.readouterr()
    assert "measuring 3 light subs" in o.err and "report only" in o.out


def test_main_no_lights(tmp_path):
    with pytest.raises(SystemExit):
        ss.main([str(tmp_path / "nothing*"), "-j", "1"])


@pytest.mark.parametrize("enc, flt", [("cp1252", "H?"), ("utf-8", "H\u03b1")])
def test_cli_redirected_encoding(tmp_path, cache, enc, flt):
    """Non-ASCII in headers must not crash output on a cp1252 console or redirected stdout."""
    d = tmp_path / "Herz \u00b0 Nebel" / NIGHT / "LIGHT"
    d.mkdir(parents=True)
    for i in range(3):
        write_xisf(d / f"s{i}.xisf", star_image(i), kw(FILTER="'H\u03b1'"))
    env = dict(os.environ, PYTHONIOENCODING=enc, SUBSELECT_CACHE=cache)
    with open(tmp_path / "out.txt", "wb") as fo, open(tmp_path / "err.txt", "wb") as fe:
        r = subprocess.run([sys.executable, SCRIPT, str(d), "-j", "2"],
                           stdout=fo, stderr=fe, env=env, timeout=300)
    out = (tmp_path / "out.txt").read_bytes().decode(enc)
    assert r.returncode == 0, (tmp_path / "err.txt").read_bytes().decode(enc, "replace")
    assert out.count(f" {flt} ") >= 3 and "report only" in out
