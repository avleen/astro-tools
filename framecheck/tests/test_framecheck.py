"""framecheck tests on synthetic data. ASTAP isn't needed: solve() is replaced by a known WCS."""
import math, os, sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import framecheck as fc

SCALE = 1.3 / 3600                  # deg/px
H, W = 400, 600


def make_wcs(ra, dec, pa):
    """TAN WCS centred on (ra, dec) with N.I.N.A.-convention position angle pa."""
    c, s = math.cos(math.radians(pa)), math.sin(math.radians(pa))
    return dict(CRPIX1=(W + 1) / 2, CRPIX2=(H + 1) / 2, CRVAL1=ra, CRVAL2=dec,
                CD1_1=-SCALE * c, CD1_2=SCALE * s, CD2_1=SCALE * s, CD2_2=SCALE * c)


def write_fits(path, **cards):
    base = [("SIMPLE", "T"), ("BITPIX", 16), ("NAXIS", 2), ("NAXIS1", W * 2), ("NAXIS2", H * 2),
            ("BZERO", 32768), ("BSCALE", 1), ("IMAGETYP", "LIGHT"), ("XPIXSZ", 3.76), ("FOCALLEN", 600)]
    lines = []
    for k, v in base + list(cards.items()):
        val = f"'{v}'" if isinstance(v, str) and v != "T" else str(v)
        lines.append(f"{k:<8}= {val:>20}".ljust(80))
    lines.append("END".ljust(80))
    head = "".join(lines)
    head += " " * (-len(head) % 2880)
    data = np.full((H * 2, W * 2), 1000 - 32768, ">i2").tobytes()
    with open(path, "wb") as f:
        f.write(head.encode() + data + b"\0" * (-len(data) % 2880))


def test_target_from_header():
    hdr = {"OBJCTRA": "'20 49 32'", "OBJCTDEC": "'+29 23 25'", "OBJCTROT": "165.0"}
    ra, dec, pa = fc.target_of(hdr)
    assert ra == pytest.approx(15 * (20 + 49 / 60 + 32 / 3600))
    assert dec == pytest.approx(29 + 23 / 60 + 25 / 3600)
    assert pa == 165.0


def test_options_override_header():
    hdr = {"OBJCTRA": "'20 49 32'", "OBJCTDEC": "'+29 23 25'", "OBJCTROT": "165.0"}
    ra, dec, pa = fc.target_of(hdr, "02:33:50", "-61:18:28", 105.0)
    assert ra == pytest.approx(15 * (2 + 33 / 60 + 50 / 3600))
    assert dec == pytest.approx(-(61 + 18 / 60 + 28 / 3600))
    assert pa == 105.0


def test_no_target():
    assert fc.target_of({}) is None
    assert fc.target_of({"OBJCTRA": "'20 49 32'", "OBJCTDEC": "'+29 23 25'"})[2] is None


@pytest.mark.parametrize("pa", [0.0, 105.0, 165.0, 285.0])
def test_position_angle_roundtrip(pa):
    assert fc.position_angle(make_wcs(10, 30, pa)) == pytest.approx(pa)


def test_on_target():
    f = fc.framing(make_wcs(312.38, 29.39, 165), (H, W), (312.38, 29.39, 165))
    assert f["offset"] == pytest.approx(0, abs=1e-6)
    assert f["dpa"] == pytest.approx(0, abs=1e-9)


def test_offset_and_shares():
    # frame centre 2' north of the target, PA 0: the target sits 2' "below" the centre
    f = fc.framing(make_wcs(100.0, 30.0 + 2 / 60, 0.0), (H, W), (100.0, 30.0, 0.0))
    assert f["offset"] == pytest.approx(2.0, rel=1e-3)
    assert abs(f["x"]) < 0.01
    assert abs(f["y"]) == pytest.approx(100 * (120 / 1.3) / H, rel=1e-3)


def test_half_turn_counts_as_same_framing():
    f = fc.framing(make_wcs(100.0, 30.0, 285.0), (H, W), (100.0, 30.0, 105.0))
    assert f["dpa"] == pytest.approx(0, abs=1e-9)
    f = fc.framing(make_wcs(100.0, 30.0, 107.5), (H, W), (100.0, 30.0, 105.0))
    assert f["dpa"] == pytest.approx(2.5)


def test_end_to_end(tmp_path, monkeypatch, capsys):
    ra, dec = 15 * (20 + 49 / 60 + 32 / 3600), 29 + 23 / 60 + 25 / 3600
    cards = dict(OBJECT="Veil Panel 1", OBJCTRA="20 49 32", OBJCTDEC="+29 23 25", OBJCTROT=165.0)
    write_fits(tmp_path / "a.fits", DATE_LOC="x", **cards)
    write_fits(tmp_path / "b.fits", **cards)
    off = {"a.fits": 0.0, "b.fits": 30 / 60}          # b is 30' off: 23% of a 400 px x 1.3" frame height

    def fake_solve(img, hdr, scale, exe, tmpdir, tag):
        assert img.shape == (H, W)                     # binned 2x2 before solving
        assert scale == pytest.approx(2 * 206.265 * 3.76 / 600)
        name = "b.fits" if hdr is b_hdr[0] else "a.fits"
        return make_wcs(ra, dec + off[name], 165.0)

    b_hdr = [None]
    real_read = fc.read_image

    def read(path):
        img, hdr = real_read(path)
        if path.endswith("b.fits"):
            b_hdr[0] = hdr
        return img, hdr

    monkeypatch.setattr(fc, "read_image", read)
    monkeypatch.setattr(fc, "solve", fake_solve)
    monkeypatch.setattr(fc, "default_astap_dir", lambda: str(tmp_path))
    monkeypatch.setattr(fc, "find_database", lambda d: "h18")
    monkeypatch.setattr(fc, "astap_exe", lambda d: "astap_cli")
    rc = fc.main([str(tmp_path), "-j", "1"])
    out = capsys.readouterr().out
    assert "Veil Panel 1: 2 subs" in out
    assert "1 flagged" in out and "off target" in out
    assert rc == 1


def test_needs_both_ra_and_dec(tmp_path):
    with pytest.raises(SystemExit):
        fc.main([str(tmp_path), "--ra", "20 49 32"])
