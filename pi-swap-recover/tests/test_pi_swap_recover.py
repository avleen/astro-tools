import os, struct, subprocess, sys, xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import pi_swap_recover as psr  # noqa: E402

SCRIPT = HERE / "pi_swap_recover.py"
NS = "{http://www.pixinsight.com/xisf}"
T0 = 1_800_000_000          # a fixed base mtime


def make_swap(d, sid, view, img, rgb=None, parts=1, mtime=T0, dirs=None, drop=()):
    """Write a synthetic PixInsight swap set for img (nc, h, w). Returns the paths."""
    img = np.asarray(img)
    nc, h, w = img.shape
    rgb = nc == 3 if rgb is None else rgb
    fl = img.dtype.kind == "f"
    hdr = bytearray(128)
    hdr[:8] = b"egamiwar"
    hdr[12], hdr[15] = int(fl), img.dtype.itemsize * 8
    struct.pack_into("<4I", hdr, 0x1C, w, h, nc, 1 if rgb else 0)
    struct.pack_into("<I", hdr, 0x68, parts)
    flat = img.astype(img.dtype.newbyteorder("<")).reshape(nc, parts, -1)    # [ch][part][slice]
    dirs = dirs or [d]
    paths = []
    for k in range(parts):
        if k in drop:
            continue
        p = Path(dirs[k % len(dirs)]) / (f"~PI~{sid}~{view}" + (f"-part{k:03d}" if k else "") + ".swp")
        p.write_bytes((bytes(hdr) if k == 0 else b"") + flat[:, k, :].tobytes())
        os.utime(p, (mtime, mtime))
        paths.append(p)
    return paths


def read_xisf(path):
    b = Path(path).read_bytes()
    assert b[:8] == b"XISF0100"
    n = struct.unpack_from("<I", b, 8)[0]
    root = ET.fromstring(b[16:16 + n].decode("utf-8"))
    im = root.find(NS + "Image")
    w, h, nc = map(int, im.get("geometry").split(":"))
    dt = {"Float32": "<f4", "Float64": "<f8", "UInt8": "u1", "UInt16": "<u2", "UInt32": "<u4"}[im.get("sampleFormat")]
    _, off, size = im.get("location").split(":")
    data = np.frombuffer(b[int(off):int(off) + int(size)], dt).reshape(nc, h, w)
    return data, im, root


def run(*args):
    return psr.main([str(a) for a in args])


def outputs(d, ext=".xisf"):
    return sorted(p for p in Path(d).iterdir() if p.suffix == ext)


@pytest.fixture
def rng():
    return np.random.default_rng(1)


def test_mono_float32_roundtrip(tmp_path, rng):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    img = rng.random((1, 40, 64), dtype=np.float32)
    make_swap(sw, "AAAAAAAAAAAAAAA1", "MasterLight_H", img, parts=8)
    assert run(sw, out) == 0
    [x] = outputs(out)
    data, im, _ = read_xisf(x)
    assert im.get("colorSpace") == "Gray" and im.get("sampleFormat") == "Float32"
    np.testing.assert_array_equal(data, img)
    assert x.name.endswith("_MasterLight_H_AAAAAAAAAAAAAAA1.xisf")


def test_rgb_float32_channel_order(tmp_path, rng):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    img = rng.random((3, 32, 48), dtype=np.float32)
    img[0] += 10; img[1] += 20; img[2] += 30           # make channels unmistakable
    make_swap(sw, "BBBBBBBBBBBBBBB2", "RGB", img, parts=8)
    assert run(sw, out) == 0
    data, im, _ = read_xisf(outputs(out)[0])
    assert im.get("colorSpace") == "RGB"
    np.testing.assert_array_equal(data, img)
    assert [int(data[c].mean()) for c in range(3)] == [10, 20, 30]


def test_uint16(tmp_path, rng):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    img = rng.integers(0, 65536, (1, 16, 24), dtype=np.uint16)
    make_swap(sw, "CCCCCCCCCCCCCCC3", "Int16", img, parts=2)
    assert run(sw, out, "--previews") == 0
    data, im, _ = read_xisf(outputs(out)[0])
    assert im.get("sampleFormat") == "UInt16" and im.get("bounds") is None
    np.testing.assert_array_equal(data, img)
    assert len(outputs(out, ".png")) == 1


def test_missing_part_and_bad_magic_skipped(tmp_path, rng, capsys):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    img = rng.random((1, 16, 16), dtype=np.float32)
    make_swap(sw, "DDDDDDDDDDDDDDD4", "Broken", img, parts=4, drop=(2,), mtime=T0)
    make_swap(sw, "DDDDDDDDDDDDDDD5", "NoMain", img, parts=4, drop=(0,), mtime=T0 + 1)
    (sw / "~PI~EEEEEEEEEEEEEEE6~Junk.swp").write_bytes(b"x" * 300)
    make_swap(sw, "FFFFFFFFFFFFFFF7", "Good", img, mtime=T0 + 2)
    assert run(sw, out) == 0
    o = capsys.readouterr().out
    assert "missing part(s) [2]" in o and "main file" in o and "not a PixInsight swap file" in o
    assert [p.name.split("_", 2)[2] for p in outputs(out)] == ["Good_FFFFFFFFFFFFFFF7.xisf"]


def test_latest_view_and_mtime_order(tmp_path, rng, capsys):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    img = rng.random((1, 8, 8), dtype=np.float32)
    # ids sort opposite to mtimes, so only mtime can give the right order
    make_swap(sw, "ZZZZZZZZZZZZZZZ1", "A", img, mtime=T0)
    make_swap(sw, "YYYYYYYYYYYYYYY2", "A", img + 1, mtime=T0 + 100)
    make_swap(sw, "XXXXXXXXXXXXXXX3", "B", img, mtime=T0 + 50)
    states = psr.find_states(str(sw))
    assert [s["id"][0] for s in states] == ["Z", "X", "Y"]
    assert run(sw, out, "--latest") == 0
    assert sorted(p.stem.split("_", 2)[2] for p in outputs(out)) == ["A_YYYYYYYYYYYYYYY2", "B_XXXXXXXXXXXXXXX3"]
    np.testing.assert_array_equal(read_xisf(next(p for p in outputs(out) if "_A_" in p.name))[0], img + 1)
    out2 = tmp_path / "out2"
    assert run(sw, out2, "--view", "A") == 0
    assert len(outputs(out2)) == 2
    out3 = tmp_path / "out3"
    assert run("--latest", sw, "--view", "A", out3) == 0          # options between positionals
    assert [p.stem.split("_", 2)[2] for p in outputs(out3)] == ["A_YYYYYYYYYYYYYYY2"]
    assert run(sw, tmp_path / "out4", "--view", "Nope") == 1


def test_previews(tmp_path, rng):
    from PIL import Image
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    make_swap(sw, "GGGGGGGGGGGGGGG1", "M", rng.random((1, 32, 40), dtype=np.float32), mtime=T0)
    make_swap(sw, "GGGGGGGGGGGGGGG2", "C", rng.random((3, 32, 40), dtype=np.float32), parts=8, mtime=T0 + 1)
    assert run(sw, out, "--previews") == 0
    pngs = outputs(out, ".png")
    assert [p.stem for p in pngs] == [p.stem for p in outputs(out)]
    assert sorted(Image.open(p).mode for p in pngs) == ["L", "RGB"]


def test_safe_name():
    assert psr.safe_name('a<b>c:d"e/f\\g|h?i*j') == "a_b_c_d_e_f_g_h_i_j"
    assert psr.safe_name("x" * 200) == "x" * 64
    assert psr.safe_name("trail. ") == "trail"
    assert psr.safe_name("") == "_"


@pytest.mark.skipif(os.name == "nt", reason="these names can't exist on Windows")
def test_invalid_view_chars_in_filename(tmp_path, rng):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    make_swap(sw, "HHHHHHHHHHHHHHH1", 'a:b?c*d"e|f<g>h\\i', rng.random((1, 8, 8), dtype=np.float32))
    assert run(sw, out) == 0
    assert outputs(out)[0].name.endswith("_a_b_c_d_e_f_g_h_i_HHHHHHHHHHHHHHH1.xisf")


def test_long_view_name_and_deep_out_dir(tmp_path, rng):
    sw = tmp_path / "swap"
    sw.mkdir()
    view = "V" * 120
    img = rng.random((1, 8, 8), dtype=np.float32)
    make_swap(sw, "IIIIIIIIIIIIIIII", view, img)
    out = tmp_path.joinpath(*["d" * 40] * 5)                       # pushes the full path past 260 chars
    assert run(sw, psr.long_path(str(out))) == 0
    [x] = [p for p in os.listdir(psr.long_path(str(out))) if p.endswith(".xisf")]
    assert x.endswith("_" + "V" * 64 + "_IIIIIIIIIIIIIIII.xisf")
    np.testing.assert_array_equal(read_xisf(psr.long_path(os.path.join(str(out), x)))[0], img)


def test_xml_escaping(tmp_path):
    p = tmp_path / "x.xisf"
    note = 'view A&B <x> "quoted" \'s C:\\Temp\\é'
    psr.write_xisf(str(p), np.zeros((1, 2, 2), np.float32), False, note)
    _, im, _ = read_xisf(p)                                     # ET.fromstring would fail on bad XML
    assert im.find(NS + "FITSKeyword").get("comment") == note


def test_ampersand_view_end_to_end(tmp_path, rng):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    make_swap(sw, "JJJJJJJJJJJJJJJ1", "R&D", rng.random((1, 8, 8), dtype=np.float32))
    assert run(sw, out) == 0
    _, im, _ = read_xisf(outputs(out)[0])
    assert "R&D" in im.find(NS + "FITSKeyword").get("comment")


def test_parts_across_several_dirs_and_windows_style_paths(tmp_path, rng):
    d1, d2, out = tmp_path / "s1", tmp_path / "s2", tmp_path / "out"
    d1.mkdir(); d2.mkdir()
    img = rng.random((3, 16, 16), dtype=np.float32)
    make_swap(d1, "KKKKKKKKKKKKKKK1", "Split", img, parts=4, dirs=[d1, d2])
    a, b = str(d1), str(d2).replace("\\", "/") + '"'               # forward slashes; cmd's stray trailing quote
    assert run(a, b, a, out) == 0                                   # duplicate dir is ignored
    np.testing.assert_array_equal(read_xisf(outputs(out)[0])[0], img)


def test_copy_to_preserves_mtime(tmp_path, rng):
    sw, safe, out = tmp_path / "swap", tmp_path / "safe", tmp_path / "out"
    sw.mkdir()
    img = rng.random((1, 16, 16), dtype=np.float32)
    make_swap(sw, "LLLLLLLLLLLLLLL1", "A", img, parts=2, mtime=T0)
    make_swap(sw, "LLLLLLLLLLLLLLL2", "A", img + 1, parts=2, mtime=T0 + 10)
    assert run(sw, out, "--copy-to", safe, "--latest") == 0
    copies = sorted(safe.iterdir())
    assert len(copies) == 4 and {int(p.stat().st_mtime) for p in copies} == {T0, T0 + 10}
    np.testing.assert_array_equal(read_xisf(outputs(out)[0])[0], img + 1)


def test_locked_file_skipped(tmp_path, rng, monkeypatch, capsys):
    import builtins
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    img = rng.random((1, 8, 8), dtype=np.float32)
    [locked] = make_swap(sw, "MMMMMMMMMMMMMMM1", "Locked", img, mtime=T0)
    make_swap(sw, "MMMMMMMMMMMMMMM2", "Fine", img, mtime=T0 + 1)
    real_open = builtins.open

    def fake_open(f, *a, **k):
        if os.path.abspath(str(f)) == str(locked):
            raise PermissionError(13, "The process cannot access the file", str(f))
        return real_open(f, *a, **k)
    monkeypatch.setattr(builtins, "open", fake_open)
    assert run(sw, out) == 0
    o = capsys.readouterr().out
    assert "locked" in o and "copy the swap files first" in o
    assert len(outputs(out)) == 1


@pytest.mark.skipif(os.name != "nt", reason="Windows byte-range lock")
def test_really_locked_file_windows(tmp_path, rng, capsys):
    import msvcrt
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    [p] = make_swap(sw, "NNNNNNNNNNNNNNN1", "Locked", rng.random((1, 8, 8), dtype=np.float32))
    with open(p, "r+b") as f:
        msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 128)
        try:
            assert run(sw, out) == 0
        finally:
            f.seek(0)
            msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 128)
    assert "SKIPPED" in capsys.readouterr().out and not outputs(out)


def test_default_swap_dirs_from_pixinsight_settings(tmp_path, monkeypatch):
    cfg = tmp_path / ("AppData" if os.name == "nt" else ".config") / "Pleiades"
    cfg.mkdir(parents=True)
    (cfg / "PixInsight.ini").write_bytes(
        b"[000]\r\nGlobal\\Preferences\\ImageWindow\\SwapDirectory001=D:/Temp\r\n"
        b"Global\\Preferences\\ImageWindow\\SwapCompression=false\r\n"
        b"Global\\Preferences\\ImageWindow\\SwapDirectory002=E:/Swap\r\n")
    monkeypatch.setenv("APPDATA", str(tmp_path / "AppData"))
    monkeypatch.setenv("HOME", str(tmp_path))
    assert psr.pi_swap_dirs() == ["D:/Temp", "E:/Swap"]
    (cfg / "PixInsight.ini").unlink()
    import tempfile
    assert psr.pi_swap_dirs() == [tempfile.gettempdir()]


def test_no_swap_files_message(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    assert run(tmp_path / "empty", tmp_path / "out") == 1
    assert "Global Preferences" in capsys.readouterr().out


def test_cp1252_redirected_stdout(tmp_path, rng):
    sw, out = tmp_path / "swap", tmp_path / "out"
    sw.mkdir()
    make_swap(sw, "OOOOOOOOOOOOOOO1", "M42_星雲", rng.random((1, 8, 8), dtype=np.float32))
    log = tmp_path / "log.txt"
    env = dict(os.environ, PYTHONIOENCODING="cp1252", PYTHONUTF8="0")
    with open(log, "wb") as f:
        r = subprocess.run([sys.executable, str(SCRIPT), str(sw), str(out), "--previews"],
                           stdout=f, stderr=subprocess.PIPE, env=env)
    assert r.returncode == 0, r.stderr.decode(errors="replace")
    text = log.read_bytes().decode("cp1252")
    assert "M42_??" in text
    assert [p.name.split("_", 2)[2] for p in outputs(out)] == ["M42_星雲_OOOOOOOOOOOOOOO1.xisf"]
