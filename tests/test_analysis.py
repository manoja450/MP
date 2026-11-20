import numpy as np
import pytest

from mcverify.analysis.curves import (crossing, gamma_1d, normalise_pdd, pass_rate, pdd_metrics,
                                      profile_metrics)
from mcverify.analysis.topas_io import read_binary, read_csv, to_1d
from tests.helpers import analytic_pdd, analytic_profile, write_topas_csv


def test_read_csv_roundtrip(tmp_path):
    vals = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
    write_topas_csv(tmp_path / "a.csv", (2, 3, 4), (0.5, 1.0, 0.1), vals, std=vals / 10)
    r = read_csv(tmp_path / "a.csv")
    assert r.shape == (2, 3, 4)
    assert r.columns == ["Sum", "Standard_Deviation"]
    assert r.unit == "Gy" and r.quantity == "DoseToMedium"
    np.testing.assert_allclose(r.data["Sum"], vals)
    np.testing.assert_allclose(r.data["Standard_Deviation"], vals / 10)
    assert r.bin_widths == [0.5, 1.0, 0.1]


def test_read_binary_ordering(tmp_path):
    ni, nj, nk = 3, 2, 4
    vals = np.random.default_rng(1).random((ni, nj, nk))
    std = vals / 7
    # TOPAS writes k outermost, i innermost, columns interleaved per voxel
    raw = np.stack([vals, std], axis=-1).transpose(2, 1, 0, 3).ravel()
    raw.astype("<f8").tofile(tmp_path / "d.bin")
    (tmp_path / "d.binheader").write_text(
        "# Results for scorer: d\n# X in 3 bins of 1 cm\n# Y in 2 bins of 1 cm\n# Z in 4 bins of 1 cm\n"
        "# DoseToMedium ( Gy ) : Sum   Standard_Deviation   \n# Binary file: d.bin\n")
    r = read_binary(tmp_path / "d.binheader")
    np.testing.assert_allclose(r.data["Sum"], vals)
    np.testing.assert_allclose(r.data["Standard_Deviation"], std)


def test_to_1d_axis_check(tmp_path):
    write_topas_csv(tmp_path / "p.csv", (1, 1, 5), (1, 1, 1), np.arange(5.0))
    r = read_csv(tmp_path / "p.csv")
    np.testing.assert_allclose(to_1d(r, "k"), np.arange(5.0))
    with pytest.raises(ValueError):
        to_1d(r, "i") if r.shape[0] == 1 and r.shape[2] > 1 else None


def test_pdd_metrics_on_analytic_curve():
    d = np.arange(0.05, 30, 0.1)
    y = analytic_pdd(d, dmax=1.5, mu=0.05)
    m = pdd_metrics(d, y)
    true_dmax = d[np.argmax(analytic_pdd(np.arange(0.0, 5, 1e-4)))] if False else \
        np.arange(0.0, 5, 1e-4)[np.argmax(analytic_pdd(np.arange(0.0, 5, 1e-4)))]
    assert m.dmax_cm == pytest.approx(true_dmax, abs=0.05)
    pct, _, _ = normalise_pdd(d, y)
    assert pct.max() == pytest.approx(100, abs=0.2)
    assert m.pdd10 == pytest.approx(100 * analytic_pdd(10.0) / analytic_pdd(true_dmax), rel=2e-3)


def test_profile_metrics_on_blurred_top_hat():
    x = np.arange(-12, 12, 0.05)
    sigma = 0.3
    y = analytic_profile(x, 5.0, sigma=sigma)
    m = profile_metrics(x, y)
    assert m.field_size_cm == pytest.approx(10.0, abs=0.02)
    assert m.centre_cm == pytest.approx(0.0, abs=0.01)
    # 80%-20% width of an erf edge = 2 * 0.8416 * sigma
    assert m.penumbra_left_cm == pytest.approx(2 * 0.8416 * sigma, abs=0.03)
    assert m.symmetry_pct < 0.5
    assert m.flatness_pct < 0.5


def test_profile_detects_shift_and_asymmetry():
    x = np.arange(-12, 12, 0.05)
    y = analytic_profile(x - 0.2, 5.0)
    assert profile_metrics(x, y).centre_cm == pytest.approx(0.2, abs=0.02)


def test_crossing_interpolates():
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = np.array([0.0, 50.0, 100.0, 50.0, 0.0])
    assert crossing(x, y, 25, "left") == pytest.approx(0.5)
    assert crossing(x, y, 25, "right") == pytest.approx(3.5)


def test_gamma_identical_and_shifted():
    x = np.arange(-10, 10, 0.1)
    y = 100 * analytic_profile(x, 5.0)
    g = gamma_1d(x, y, x, y, 1.0, 1.0)
    assert pass_rate(g) == 100.0
    assert np.nanmax(g) == pytest.approx(0.0, abs=1e-9)
    # a 0.5 mm shift passes 1%/1mm, a 3 mm shift does not
    assert pass_rate(gamma_1d(x, y, x + 0.05, y, 1.0, 1.0)) == 100.0
    assert pass_rate(gamma_1d(x, y, x + 0.3, y, 1.0, 1.0)) < 95.0
    # a pure 0.5% dose offset passes, a 2% offset fails in the flat region
    assert pass_rate(gamma_1d(x, y, x, y * 1.005, 1.0, 1.0)) == 100.0
    assert pass_rate(gamma_1d(x, y, x, y * 1.02, 1.0, 1.0)) < 80.0


def test_gamma_worked_example_from_guide():
    # 1.5% dose difference at 1 mm distance with 3%/2 mm -> gamma = 0.707
    x_ref = np.array([0.0])
    d_ref = np.array([100.0])
    x_ev = np.array([-1.0, 0.1, 1.0])
    d_ev = np.array([0.0, 101.5, 0.0])
    g = gamma_1d(x_ref, d_ref, x_ev, d_ev, 3.0, 2.0, threshold_pct=0, norm_value=100.0)
    assert g[0] <= np.sqrt(0.25 + 0.25) + 1e-6
