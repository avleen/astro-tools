"""subdepth tests on synthetic data: star database tiles, WCS, sky noise, zero point, end to end.
ASTAP itself is not needed: the plate solve is replaced by the known WCS."""
import csv, os

import numpy as np
import pytest

import subdepth as sd

SCALE = 1.3                         # arcsec/px
RA0, DEC0 = 40.0, 61.0              # field centre
H, W = 512, 640


# ---------------------------------------------------------------- helpers
def encode_tile(path, ra, dec, mag):
    """Write a .1476 file the way ASTAP stores H17/H18 (5-byte records, magnitude groups)."""
    head = b"synthetic test tile".ljust(109) + bytes([5])
    recs = []
    dec_raw = np.round(np.asarray(dec) / 90 * (2 ** 23 - 1)).astype(np.int64)
    ra_raw = np.round(np.asarray(ra) / 360 * (2 ** 24 - 1)).astype(np.int64) % (2 ** 24 - 1)
    m10 = np.round(np.asarray(mag) * 10).astype(np.int64)
    order = np.lexsort((dec_raw >> 16, m10))
    last = None
    for i in order:
        key = (m10[i], dec_raw[i] >> 16)
        if key != last:
            recs.append(bytes([255, 255, 255, int(key[1] + 128), int(key[0] + 16)]))
            last = key
        r, d = ra_raw[i], dec_raw[i]
        recs.append(bytes([r & 255, (r >> 8) & 255, (r >> 16) & 255, d & 255, (d >> 8) & 255]))
    with open(path, "wb") as f:
        f.write(head + b"".join(recs))


def tile_name(ra, dec):
    for r in range(1, 37):
        lo, hi = sd.ring_dec(r)
        if lo <= dec < hi:
            n = sd.NCELL[r - 1]
            return f"h18_{r:02d}{int(ra // (360 / n)) + 1:02d}.1476"


def make_wcs(rot_deg=20.0):
    c, s = np.cos(np.radians(rot_deg)), np.sin(np.radians(rot_deg))
    k = SCALE / 3600
    return dict(CRPIX1=W / 2 + 0.5, CRPIX2=H / 2 + 0.5, CRVAL1=RA0, CRVAL2=DEC0,
                CD1_1=-k * c, CD1_2=k * s, CD2_1=k * s, CD2_2=k * c)


def pixel_to_world(w, x, y):
    """Inverse TAN, for building the synthetic catalogue."""
    cd = np.array([[w["CD1_1"], w["CD1_2"]], [w["CD2_1"], w["CD2_2"]]])
    u, v = x - (w["CRPIX1"] - 1), y - (w["CRPIX2"] - 1)
    xi, eta = np.radians(cd[0, 0] * u + cd[0, 1] * v), np.radians(cd[1, 0] * u + cd[1, 1] * v)
    a0, d0 = np.radians(w["CRVAL1"]), np.radians(w["CRVAL2"])
    den = np.cos(d0) - eta * np.sin(d0)
    ra = a0 + np.arctan2(xi, den)
    dec = np.arctan2(np.sin(d0) + eta * np.cos(d0), np.hypot(xi, den))
    return np.degrees(ra) % 360, np.degrees(dec)


def write_fits(path, img, **cards):
    """Uncompressed 16-bit FITS light frame with header cards."""
    base = [("SIMPLE", "T"), ("BITPIX", 16), ("NAXIS", 2), ("NAXIS1", img.shape[1]), ("NAXIS2", img.shape[0]),
            ("BZERO", 32768), ("BSCALE", 1)]
    lines = []
    for k, v in base + list(cards.items()):
        val = f"'{v}'" if isinstance(v, str) and v != "T" else str(v)
        lines.append(f"{k:<8}= {val:>20}".ljust(80))
    lines.append("END".ljust(80))
    head = "".join(lines)
    head += " " * (-len(head) % 2880)
    data = (np.clip(np.round(img), 0, 65535) - 32768).astype(">i2").tobytes()
    with open(path, "wb") as f:
        f.write(head.encode() + data + b"\0" * (-len(data) % 2880))


def synth_field(rng, zp_e, sky_e, gain, exp, rn=1.5, bias=500.0, n=150, fwhm=3.0):
    """Light frame in ADU plus its catalogue (ra, dec, mag) and WCS."""
    w = make_wcs()
    x = rng.uniform(30, W - 30, n)
    y = rng.uniform(30, H - 30, n)
    mag = rng.uniform(10.5, 14.5, n)
    ra, dec = pixel_to_world(w, x, y)
    yy, xx = np.mgrid[0:H, 0:W]
    model = np.full((H, W), sky_e, float)
    sig = fwhm / 2.3548
    e_tot = 10 ** (-0.4 * (mag - zp_e)) * exp
    for xi, yi, e in zip(x, y, e_tot):
        sl = (slice(int(yi) - 12, int(yi) + 13), slice(int(xi) - 12, int(xi) + 13))
        model[sl] += e / (2 * np.pi * sig ** 2) * np.exp(-((xx[sl] - xi) ** 2 + (yy[sl] - yi) ** 2) / (2 * sig ** 2))
    e = rng.poisson(model) + rng.normal(0, rn, model.shape)
    return e / gain + bias, (ra, dec, mag), w


# ---------------------------------------------------------------- star database
def test_tile_roundtrip(tmp_path):
    rng = np.random.default_rng(1)
    ra = rng.uniform(0, 5.2, 500)
    dec = rng.uniform(-5.14, 5.0, 500)          # crosses dec 0: negative high bytes
    mag = rng.uniform(4, 18, 500)
    p = tmp_path / "t.1476"
    encode_tile(p, ra, dec, mag)
    r2, d2, m2 = sd.read_tile(str(p))
    o1, o2 = np.lexsort((ra, np.round(mag, 1))), np.lexsort((r2, m2))
    assert len(r2) == 500
    assert np.allclose(np.sort(ra), np.sort(r2), atol=1e-4)
    assert np.allclose(np.sort(dec), np.sort(d2), atol=1e-4)
    assert np.allclose(np.round(mag[o1], 1), m2[o2], atol=1e-6)


def test_read_tile_rejects_other_formats(tmp_path):
    p = tmp_path / "x.1476"
    p.write_bytes(b"x" * 109 + bytes([11]) + b"\0" * 22)
    with pytest.raises(ValueError):
        sd.read_tile(str(p))


def test_rings_cover_the_sky():
    assert sum(sd.NCELL) == 1476
    assert sd.ring_dec(1)[0] == -90 and sd.ring_dec(36)[1] == 90
    for r in range(1, 36):
        assert sd.ring_dec(r)[1] == pytest.approx(sd.ring_dec(r + 1)[0])
    assert sd.ring_dec(18)[1] == pytest.approx(0, abs=1e-9)


def test_catalog_reads_only_nearby_cells(tmp_path):
    near = tile_name(RA0, DEC0)
    far = tile_name(RA0 + 90, DEC0)
    encode_tile(tmp_path / near, [RA0, RA0 + 0.5], [DEC0, DEC0 + 0.2], [12.0, 13.0])
    encode_tile(tmp_path / far, [RA0 + 90], [DEC0], [11.0])
    ra, dec, mag = sd.catalog(str(tmp_path), "h18", RA0, DEC0, 1.0)
    assert sorted(mag.tolist()) == [12.0, 13.0]


# ---------------------------------------------------------------- WCS
def test_wcs_roundtrip():
    w = make_wcs(35)
    x = np.array([0.0, 100.0, 400.5, W - 1.0])
    y = np.array([0.0, 300.0, 17.25, H - 1.0])
    ra, dec = pixel_to_world(w, x, y)
    x2, y2, front = sd.world_to_pixel(w, ra, dec)
    assert front.all()
    assert np.allclose(x, x2, atol=1e-6) and np.allclose(y, y2, atol=1e-6)
    assert sd.wcs_scale(w) == pytest.approx(SCALE)


# ---------------------------------------------------------------- sky and electrons
@pytest.mark.parametrize("sky_e", [8.0, 100.0])
def test_sky_from_dark_and_noise_agree(sky_e):
    rng = np.random.default_rng(2)
    gain, rn, bias = 0.25, 1.5, 500.0
    img = (rng.poisson(sky_e, (1024, 1024)) + rng.normal(0, rn, (1024, 1024))) / gain + bias
    sky_adu, sig_adu = sd.sky_stats(img.astype(np.float32))
    m = dict(bayer=False, sky_adu=sky_adu, sig_adu=sig_adu, zp_adu=None, scale=1.0)
    with_dark = sd.derive(m, 60, gain, rn, bias)
    from_noise = sd.derive(m, 60, gain, rn, None)
    assert with_dark["sky_method"] == "dark" and from_noise["sky_method"] == "noise"
    assert with_dark["sky_e"] == pytest.approx(sky_e, rel=0.05, abs=0.5)
    assert from_noise["sky_e"] == pytest.approx(sky_e, rel=0.08, abs=0.8)
    assert with_dark["gain_chk"] == pytest.approx(gain, rel=0.1)


def test_superpixel_sums_bayer_cells():
    a = np.arange(16, dtype=float).reshape(4, 4)
    assert sd.superpixel(a).tolist() == [[0 + 1 + 4 + 5, 2 + 3 + 6 + 7], [8 + 9 + 12 + 13, 10 + 11 + 14 + 15]]


def test_mu_from_zero_point_and_sky():
    m = dict(bayer=False, sky_adu=600.0, sig_adu=1.0, zp_adu=20.0, scale=2.0)
    o = sd.derive(m, 100, 0.5, 1.0, 500.0)       # 100 ADU * 0.5 = 50 e- in 100 s = 0.5 e-/s/px = 0.125 /arcsec^2
    assert o["zp"] == pytest.approx(20 + 2.5 * np.log10(0.5))
    assert o["rate"] == pytest.approx(0.125)
    assert o["mu"] == pytest.approx(o["zp"] - 2.5 * np.log10(0.125))


# ---------------------------------------------------------------- effective time
def test_effective_time():
    assert sd.effective_time(300, 1.0, 0.2, 0.2, 1.0, 0.0) == pytest.approx(300)
    assert sd.effective_time(300, 0.5, 0.2, 0.2, 1.0, 0.0) == pytest.approx(75)       # T^2
    assert sd.effective_time(300, 1.0, 0.8, 0.2, 1.0, 0.0) == pytest.approx(75)       # 4x the sky
    with_rn = sd.effective_time(300, 1.0, 0.2, 0.2, 1.0, 1.25)
    assert with_rn == pytest.approx(300 * 60 / (60 + 1.5625))                           # sky 60 e- vs RN^2


def test_add_depth_references_per_filter():
    subs = [dict(filter="H", exp=300, zp=15.0 - 0.01 * i, rate=0.1 * (1 + i), scale=1.0, bayer=False) for i in range(20)]
    subs += [dict(filter="O", exp=300, zp=14.0, rate=0.05, scale=1.0, bayer=False)]
    sd.add_depth(subs, 0.0)
    h = subs[:20]
    assert all(0 < s["T"] <= 1 for s in h)
    assert h[0]["t_eff"] == pytest.approx(300, rel=0.05)          # clearest, darkest sub: nearly its full time
    assert h[-1]["t_eff"] < h[0]["t_eff"] / 5                       # 20x the sky
    assert subs[-1]["T"] == 1 and subs[-1]["t_eff"] == pytest.approx(300)   # O has its own reference


# ---------------------------------------------------------------- zero point and end to end
def test_zero_point_recovered(tmp_path):
    rng = np.random.default_rng(3)
    zp_e, gain, exp = 18.0, 0.25, 60.0
    img, (ra, dec, mag), w = synth_field(rng, zp_e, 20.0, gain, exp)
    encode_tile(tmp_path / tile_name(RA0, DEC0), ra, dec, mag)
    stars, _ = sd.find_stars(img.astype(np.float32))
    zp_adu, scatter, n = sd.zero_point(img.astype(np.float32), stars, w, str(tmp_path), "h18", exp)
    assert n > 50
    # aperture of 1.5 FWHM holds ~99.8% of a Gaussian; catalogue rounds to 0.1 mag
    assert zp_adu + 2.5 * np.log10(gain) == pytest.approx(zp_e, abs=0.03)
    assert scatter < 0.05


def test_end_to_end(tmp_path, monkeypatch):
    rng = np.random.default_rng(4)
    astap = tmp_path / "astap"
    astap.mkdir()
    (astap / "astap_cli").write_text("stand-in")
    (astap / "h18_0101.1476").write_bytes(b"x" * 109 + bytes([5]))
    lights = tmp_path / "lights" / "NIGHT_2026-10-22"
    lights.mkdir(parents=True)
    cat = None
    # three subs: clear + dark, hazy (T = 0.5), bright sky (4x)
    for i, (zp_e, sky) in enumerate([(18.0, 20.0), (18.0 - 2.5 * np.log10(2), 20.0), (18.0, 80.0)]):
        img, cat_i, w = synth_field(np.random.default_rng(10), zp_e, sky, 0.25, 60.0)
        img += rng.normal(0, 0.1, img.shape)
        cat = cat or cat_i
        write_fits(str(lights / f"L{i}.fits"), img, IMAGETYP="Light Frame", EXPTIME=60.0, FILTER="L", EGAIN=0.25,
                   XPIXSZ=3.76, FOCALLEN=round(206.265 * 3.76 / SCALE, 3), RA=RA0, DEC=DEC0,
                   **{"DATE-OBS": f"2026-10-23T03:0{i}:00"})
    encode_tile(astap / tile_name(RA0, DEC0), *cat)
    monkeypatch.setattr(sd, "solve", lambda *a, **k: make_wcs())
    out = tmp_path / "depth.csv"
    sd.main([str(tmp_path / "lights"), "--astap", str(astap), "--read-noise", "1.5", "-j", "1", "--csv", str(out)])
    rows = {r["file"]: r for r in csv.DictReader(open(out, encoding="utf-8-sig"))}
    assert all(r["solved"] == "True" and r["sky_method"] == "noise" for r in rows.values())
    zp = [float(rows[f"L{i}.fits"]["zp"]) for i in range(3)]
    assert zp[0] == pytest.approx(18.0, abs=0.05) and zp[1] == pytest.approx(18.0 - 0.753, abs=0.05)
    T = [float(rows[f"L{i}.fits"]["T"]) for i in range(3)]
    assert T[1] == pytest.approx(0.5, abs=0.04)
    sky = [float(rows[f"L{i}.fits"]["sky_e"]) for i in range(3)]
    assert sky[2] / sky[0] == pytest.approx(4, rel=0.1)
    te = [float(rows[f"L{i}.fits"]["t_eff"]) for i in range(3)]
    assert te[1] / te[0] == pytest.approx(0.25, rel=0.15)
    assert te[2] / te[0] == pytest.approx(0.25, rel=0.15)
    assert not list((tmp_path / "lights").rglob("*.ini")) and not list((tmp_path / "lights").rglob("*.wcs"))


def test_gain_per_setting():
    assert sd.gain_for(["100=0.25", "0=0.79"], {"GAIN": "100."}) == 0.25
    assert sd.gain_for(["100=0.25", "0=0.79"], {"GAIN": "0"}) == 0.79
    assert sd.gain_for(["100=0.25"], {"GAIN": "120"}) is None
    assert sd.gain_for(["0.5", "100=0.25"], {"GAIN": "120"}) == 0.5
    assert sd.gain_for([], {}) is None


def test_adu_step():
    rng = np.random.default_rng(5)
    raw14 = rng.integers(0, 16384, (64, 64)).astype(np.float32)
    assert sd.adu_step(raw14 * 4) == 4
    assert sd.adu_step(raw14 * 2 + 1) == 1
    assert sd.adu_step(raw14 * 4 + 1) == 1
    assert sd.adu_step(raw14 / 7) == 1                 # float data


def test_egain_divided_for_scaled_14bit(tmp_path, monkeypatch):
    """A 14-bit frame stored x4 with EGAIN per native ADU gives the same electrons as the true file gain."""
    astap = tmp_path / "astap"
    astap.mkdir()
    (astap / "astap_cli").write_text("stand-in")
    (astap / "h18_0101.1476").write_bytes(b"x" * 109 + bytes([5]))
    img, cat, w = synth_field(np.random.default_rng(6), 18.0, 50.0, 1.0, 60.0, bias=200)
    img = np.round(img) * 4                                  # native ADU at 1 e-/ADU, stored x4
    encode_tile(astap / tile_name(RA0, DEC0), *cat)
    d = tmp_path / "lights"
    d.mkdir()
    write_fits(str(d / "a.fits"), img, IMAGETYP="LIGHT", EXPTIME=60.0, EGAIN=1.0, XPIXSZ=3.76,
               FOCALLEN=round(206.265 * 3.76 / SCALE, 3), RA=RA0, DEC=DEC0)
    monkeypatch.setattr(sd, "solve", lambda *a, **k: make_wcs())
    out = tmp_path / "o.csv"
    sd.main([str(d), "--astap", str(astap), "--read-noise", "1.5", "-j", "1", "--csv", str(out)])
    r = next(csv.DictReader(open(out, encoding="utf-8-sig")))
    assert float(r["zp"]) == pytest.approx(18.0, abs=0.05)
    assert float(r["sky_e"]) == pytest.approx(50.0, rel=0.1)


# ---------------------------------------------------------------- bias level without darks
def test_estimated_pedestal_from_offset():
    assert sd.estimated_pedestal({"INSTRUME": "ZWO ASI2600MM Pro", "OFFSET": "50"}, 1) == pytest.approx(500)
    assert sd.estimated_pedestal({"INSTRUME": "ZWO CCD ASI533MC Pro", "OFFSET": "70."}, 4) == pytest.approx(2744)
    assert sd.estimated_pedestal({"INSTRUME": "QHY268M", "OFFSET": "50"}, 1) is None
    assert sd.estimated_pedestal({"INSTRUME": "ZWO ASI2600MM Pro"}, 1) is None


def test_estimated_bias_used_only_when_noise_agrees():
    rng = np.random.default_rng(7)
    gain, rn, bias, sky_e = 0.25, 1.5, 500.0, 40.0
    img = (rng.poisson(sky_e, (1024, 1024)) + rng.normal(0, rn, (1024, 1024))) / gain + bias
    sky_adu, sig_adu = sd.sky_stats(img.astype(np.float32))
    m = dict(bayer=False, sky_adu=sky_adu, sig_adu=sig_adu, zp_adu=None, scale=1.0)

    def sub(night, ped):
        s = dict(night=night, pedestal=ped, bias_source="estimated bias")
        s.update(sd.derive(m, 60, gain, rn, ped, "estimated bias"))
        s["_noise"] = sd.derive(m, 60, gain, rn, None)
        return s
    subs = [sub("good", bias), sub("good", bias), sub("wrong", bias - 40), sub("wrong", bias - 40)]  # 10 e- off
    fell = sd.check_estimates(subs)
    assert list(fell) == [("wrong", bias - 40)] and fell[("wrong", bias - 40)] > sd.AGREE
    good, wrong = subs[0], subs[2]
    assert good["sky_method"] == "estimated bias" and good["bias_noise"] == pytest.approx(1, abs=0.08)
    assert good["sky_e"] == pytest.approx(sky_e, rel=0.05)
    assert wrong["sky_method"] == "noise" and wrong["sky_e"] == pytest.approx(sky_e, rel=0.08)
    assert wrong["bias_noise"] > sd.AGREE and wrong["pedestal"] == bias - 40      # kept for the record
    dark = sd.derive(m, 60, gain, rn, bias - 40, "dark")                           # a real dark is used, but flagged
    assert dark["sky_method"] == "dark" and not sd.agrees(dark["bias_noise"])


def _one_light(tmp_path, name, **cards):
    astap = tmp_path / "astap"
    if not astap.exists():
        astap.mkdir()
        (astap / "astap_cli").write_text("stand-in")
        (astap / "h18_0101.1476").write_bytes(b"x" * 109 + bytes([5]))
    img, cat, w = synth_field(np.random.default_rng(8), 18.0, 30.0, 0.25, 60.0, rn=1.25, bias=500.0)
    encode_tile(astap / tile_name(RA0, DEC0), *cat)
    d = tmp_path / name
    d.mkdir()
    write_fits(str(d / "a.fits"), img, IMAGETYP="LIGHT", EXPTIME=60.0, EGAIN=0.25, XPIXSZ=3.76,
               FOCALLEN=round(206.265 * 3.76 / SCALE, 3), RA=RA0, DEC=DEC0, **cards)
    return astap, d


@pytest.mark.parametrize("offset, method", [(50, "estimated bias"), (60, "noise")])
def test_bias_from_offset_end_to_end(tmp_path, monkeypatch, capsys, offset, method):
    astap, d = _one_light(tmp_path, "lights", INSTRUME="ZWO ASI2600MM Pro", OFFSET=offset)
    monkeypatch.setattr(sd, "solve", lambda *a, **k: make_wcs())
    out = tmp_path / "o.csv"
    sd.main([str(d), "--astap", str(astap), "--read-noise", "1.25", "-j", "1", "--csv", str(out)])
    r = next(csv.DictReader(open(out, encoding="utf-8-sig")))
    assert r["sky_method"] == method
    assert float(r["sky_e"]) == pytest.approx(30.0, rel=0.08)
    assert float(r["pedestal"]) == pytest.approx(10 * offset)
    warned = "disagree" in capsys.readouterr().out
    assert warned == (method == "noise")


def test_pedestal_option(tmp_path, monkeypatch):
    astap, d = _one_light(tmp_path, "lights")
    monkeypatch.setattr(sd, "solve", lambda *a, **k: make_wcs())
    out = tmp_path / "o.csv"
    sd.main([str(d), "--astap", str(astap), "--read-noise", "1.25", "--pedestal", "500", "-j", "1", "--csv", str(out)])
    r = next(csv.DictReader(open(out, encoding="utf-8-sig")))
    assert r["sky_method"] == "pedestal" and float(r["sky_e"]) == pytest.approx(30.0, rel=0.05)
    assert float(r["gain_chk"]) == pytest.approx(0.25, rel=0.1)


@pytest.mark.parametrize("v, hours, want", [
    ("20 49 32", True, 15 * (20 + 49 / 60 + 32 / 3600)),
    ("+29 23 25", False, 29 + 23 / 60 + 25 / 3600),
    ("-61 18 28", False, -(61 + 18 / 60 + 28 / 3600)),
    ("-0 30 00", False, -0.5),
    ("- 5 30 00", False, -5.5),
])
def test_hms_signs(v, hours, want):
    assert sd._hms(v, hours) == pytest.approx(want)
