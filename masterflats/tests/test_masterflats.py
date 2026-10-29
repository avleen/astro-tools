import os
import re
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "masterflats.py")
sys.path.insert(0, os.path.dirname(HERE))
import masterflats as mf  # noqa: E402
from subselect import _xisf_xml, read_header, read_image  # noqa: E402  (masterflats put subselect on the path)

NIGHT = "NIGHT_2026-10-22"
OFFSET = 800.0
N = 64
FK = {"IMAGETYP": "'FLAT'", "FILTER": "'H'", "EXPTIME": "10.0", "SET-TEMP": "-10.0", "GAIN": "100",
      "DATE-LOC": "'2026-10-23T06:35:48.000'", "ROTATANG": "225.41"}


def write_xisf(path, img, kw):
    data = img.astype("<u2").tobytes()
    keys = "".join(f'<FITSKeyword name="{k}" value="{v}" comment=""/>' for k, v in kw.items())

    def xml(off):
        return (f'<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0"><Image geometry="{N}:{N}:1" '
                f'sampleFormat="UInt16" colorSpace="Gray" location="attachment:{off:010d}:{len(data)}">'
                f'{keys}</Image></xisf>').encode()
    off = 16 + len(xml(0))
    with open(path, "wb") as f:
        f.write(b"XISF0100" + len(xml(off)).to_bytes(4, "little") + b"\0" * 4 + xml(off) + data)


def vignette():
    yy, xx = np.mgrid[:N, :N]
    return 1 - 0.4 * (((xx - N / 2) ** 2 + (yy - N / 2) ** 2) / (N / 2) ** 2)     # 1.0 centre, 0.2 corners


def flat_frame(seed, level=25000.0):
    rng = np.random.default_rng(seed)
    img = OFFSET + level * vignette()
    return np.clip(rng.normal(img, np.sqrt(img)), 0, 65535).astype(np.uint16)


def add_flats(d, flt, n, start=0, **kw):
    d.mkdir(parents=True, exist_ok=True)
    for i in range(start, start + n):
        write_xisf(d / f"flat_{flt}_{i:04d}.xisf", flat_frame(100 * ord(flt[0]) + i, 20000 + 500 * i),
                   dict(FK, FILTER=f"'{flt}'", **kw))


def assert_vignette(img, tol=0.05):
    v = vignette()
    rel = (img / img.mean()) / (v / v.mean())
    assert np.abs(rel - 1).max() < tol and np.abs(rel - 1).mean() < 0.01


@pytest.fixture
def rig(tmp_path):
    """src/NIGHT/FLAT/{H,O} with 5 flats each, and a master-dark library with 60 s and 180 s darks."""
    src = tmp_path / "src"
    add_flats(src / NIGHT / "FLAT" / "H", "H", 5)
    add_flats(src / NIGHT / "FLAT" / "O", "O", 5)
    darks = tmp_path / "Darks - -10C, Gain 100" / "master"
    darks.mkdir(parents=True)
    for e in (60, 180):
        write_xisf(darks / f"{e}s.xisf", np.full((N, N), OFFSET + 2, np.uint16),
                   {"IMAGETYP": "'Master Dark'", "EXPTIME": f"{e}.0", "SET-TEMP": "-10.0"})
    return src, str(darks)


def master(src, flt, night=NIGHT):
    return src / night / "FLAT" / flt / f"MasterFlat_{night[6:]}_{flt}_10s.xisf"


def run(capsys, *args):
    rc = mf.main([str(a) for a in args])
    return rc, capsys.readouterr().out


def test_integrate_shape_and_rejection(tmp_path):
    paths = []
    for i in range(7):
        img = flat_frame(i, 15000 + 1000 * i)
        if i == 3:
            img[10, 20] = 65535                         # cosmic ray / hot pixel in one frame
        write_xisf(tmp_path / f"f{i}.xisf", img, FK)
        paths.append(str(tmp_path / f"f{i}.xisf"))
    m = mf.integrate(paths, np.full((N, N), OFFSET, np.float32), jobs=2)
    assert_vignette(m)                                  # vignetting reproduced, offset removed
    v = vignette()
    assert abs((m[10, 20] / m.mean()) / (v[10, 20] / v.mean()) - 1) < 0.03     # outlier clipped
    assert abs(m.mean() / (np.mean([15000 + 1000 * i for i in range(7)]) * v.mean()) - 1) < 0.01


def test_build_then_idempotent_then_new_flat(rig, capsys):
    src, darks = rig
    rc, o = run(capsys, src, "--darks", darks)
    assert rc == 0 and o.count(" built MasterFlat_") == 2 and "60s.xisf (nearest exposure)" in o
    for flt in "HO":
        img, h = read_image(str(master(src, flt)))
        assert h["IMAGETYP"] == "Master Flat" and h["FILTER"] == flt and h["NCOMBINE"] == "5"
        assert h["ROTATANG"] == "225.41" and float(h["EXPTIME"]) == 10 and h["GAIN"] == "100"
        assert_vignette(img)
    m0 = os.stat(master(src, "H")).st_mtime_ns

    rc, o = run(capsys, src, "--darks", darks)
    assert o.count("up to date:") == 2 and " built " not in o
    assert os.stat(master(src, "H")).st_mtime_ns == m0

    add_flats(src / NIGHT / "FLAT" / "H", "H", 1, start=5)             # a new flat rebuilds that master only
    rc, o = run(capsys, src, "--darks", darks)
    assert o.count(" built MasterFlat_") == 1 and o.count("up to date:") == 1
    assert read_header(str(master(src, "H")))["NCOMBINE"] == "6"


def test_dry_run_changes_nothing(rig, capsys):
    src, darks = rig
    rc, o = run(capsys, src, "--darks", darks, "-n", "--delete-raw")
    assert o.count("would build") == 2 and "2 to build" in o
    assert not master(src, "H").exists() and len(os.listdir(src / NIGHT / "FLAT" / "H")) == 5


def test_flat_darks_preferred_and_new_dark_rebuilds(rig, capsys):
    src, darks = rig
    run(capsys, src, "--darks", darks)
    fd = src / NIGHT / "DARK" / "H"
    fd.mkdir(parents=True)
    for i in range(3):
        write_xisf(fd / f"fd{i}.xisf", np.full((N, N), OFFSET, np.uint16), dict(FK, IMAGETYP="'DARK'"))
    rc, o = run(capsys, src, "--darks", darks)
    assert o.count(" built MasterFlat_") == 2 and o.count("3 flat-dark(s) from the same night") == 2


def test_exact_dark_and_bias_preference(rig, capsys):
    src, darks = rig
    write_xisf(os.path.join(darks, "bias.xisf"), np.full((N, N), OFFSET, np.uint16),
               {"IMAGETYP": "'Master Bias'", "EXPTIME": "0.0", "SET-TEMP": "-10.0"})
    assert "master bias bias.xisf" in run(capsys, src, "--darks", darks)[1]
    write_xisf(os.path.join(darks, "10s.xisf"), np.full((N, N), OFFSET + 1, np.uint16),
               {"IMAGETYP": "'Master Dark'", "EXPTIME": "10.0", "SET-TEMP": "-10.0"})
    o = run(capsys, src, "--darks", darks)[1]
    assert o.count(" built MasterFlat_") == 2 and "master dark 10s.xisf" in o


def test_wrong_temperature_dark_not_used(rig, tmp_path, capsys):
    src, _ = rig
    warm = tmp_path / "warm"
    warm.mkdir()
    write_xisf(warm / "10s.xisf", np.full((N, N), OFFSET, np.uint16),
               {"IMAGETYP": "'Master Dark'", "EXPTIME": "10.0", "SET-TEMP": "0.0"})
    rc, o = run(capsys, src, "--darks", warm)
    assert o.count("NOT BUILT: no flat-dark") == 2 and not master(src, "H").exists()


def test_delete_raw_keeps_master(rig, capsys):
    src, darks = rig
    rc, o = run(capsys, src, "--darks", darks, "--delete-raw")
    assert o.count("deleted 5 raw flats") == 2
    assert os.listdir(src / NIGHT / "FLAT" / "H") == [master(src, "H").name]
    rc, o = run(capsys, src, "--darks", darks)                          # nothing raw left: nothing to do
    assert "nothing to do" in o and master(src, "H").exists()


def test_partial_delete_keeps_master(rig, capsys):
    src, darks = rig
    run(capsys, src, "--darks", darks)
    m0 = os.stat(master(src, "H")).st_mtime_ns
    for f in sorted(x for x in os.listdir(src / NIGHT / "FLAT" / "H") if x.startswith("flat_"))[:3]:
        os.remove(src / NIGHT / "FLAT" / "H" / f)                      # 2 raw + master left
    rc, o = run(capsys, src, "--darks", darks)
    assert o.count("up to date:") == 2 and os.stat(master(src, "H")).st_mtime_ns == m0


def test_groups_split_and_skips(rig, capsys):
    src, darks = rig
    add_flats(src / NIGHT / "FLAT" / "H", "H", 3, start=10, ROTATANG="77.24")    # same folder, other rotation
    add_flats(src / "NIGHT_2026-10-25" / "FLAT" / "R", "R", 2)                  # too few
    add_flats(src / "_selected" / NIGHT / "FLAT", "S", 5)                       # subselect's copies: ignored
    rc, o = run(capsys, src, "--darks", darks)
    assert o.count(" built MasterFlat_") == 3 and "fewer than 3 flats" in o and "   S " not in o
    names = sorted(os.listdir(src / NIGHT / "FLAT" / "H"))
    assert f"MasterFlat_{NIGHT[6:]}_H_10s.xisf" in names and f"MasterFlat_{NIGHT[6:]}_H_10s_2.xisf" in names


def test_cli_wildcard_and_cp1252(rig, tmp_path):
    src, darks = rig
    (src / NIGHT / "FLAT" / "H" / "x.xisf").write_bytes(b"junk")     # unreadable file: skipped, no crash
    env = dict(os.environ, PYTHONIOENCODING="cp1252")
    with open(tmp_path / "log.txt", "w") as out:
        r = subprocess.run([sys.executable, SCRIPT, str(src / "NIGHT_*"), "--darks", darks],
                           stdout=out, stderr=subprocess.PIPE, env=env, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "log.txt").read_text().count(" built MasterFlat_") == 2 and "skipping" in r.stderr


def test_foreign_master_not_overwritten(rig, capsys):
    src, darks = rig
    write_xisf(master(src, "H"), np.full((N, N), 30000, np.uint16), dict(FK, IMAGETYP="'Master Flat'"))
    before = master(src, "H").read_bytes()
    rc, o = run(capsys, src, "--darks", darks)
    assert "wasn't made by masterflats" in o and master(src, "H").read_bytes() == before


def image_tag(path):
    with open(path, "rb") as f:
        return re.search(r"<Image\b[^>]*>", _xisf_xml(f)).group(0)


def test_master_is_zstd_shuffled(rig, capsys):
    src, darks = rig
    run(capsys, src, "--darks", darks)
    tag = image_tag(master(src, "H"))
    assert f'compression="zstd+sh:{N * N * 4}:4"' in tag
    assert int(re.search(r'attachment:\d+:(\d+)', tag).group(1)) < N * N * 4


def test_old_uncompressed_master_compressed_not_rebuilt(rig, capsys):
    src, darks = rig
    run(capsys, src, "--darks", darks)
    path = master(src, "H")
    img0, _ = read_image(str(path))
    with open(path, "rb") as f:                         # rewrite it as an older, uncompressed version would have
        xml = _xisf_xml(f)
    tag = image_tag(path)
    data = np.ascontiguousarray(img0 / 65535, "<f4").tobytes()
    plain = re.sub(r' compression="[^"]*"', "", tag)
    mf._write_xisf(path, lambda o: xml.replace(tag, re.sub(r'location="attachment:\d+:\d+"',
                                                           f'location="attachment:{o}:{len(data)}"', plain)).encode(), data)
    assert "compression=" not in image_tag(path)
    for f in sorted(x for x in os.listdir(path.parent) if x.startswith("flat_"))[:2]:
        os.remove(path.parent / f)                      # raw flats partly gone: a rebuild would lose frames

    o = run(capsys, src, "--darks", darks, "-n")[1]
    assert "would compress it" in o and "compression=" not in image_tag(path)
    o = run(capsys, src, "--darks", darks)[1]
    assert o.count("up to date:") == 2 and "compressed" in o and " built " not in o
    img1, h = read_image(str(path))
    assert "compression=" in image_tag(path) and h["NCOMBINE"] == "5" and np.array_equal(img0, img1)
    assert "compressed" not in run(capsys, src, "--darks", darks)[1]     # and only once
