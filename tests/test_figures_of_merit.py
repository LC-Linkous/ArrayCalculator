#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_figures_of_merit.py
#
#   Tests for the shared figure-of-merit helpers in array_common.py
#   (directivity_integral, hpbw_deg, peak_sidelobe) and for the plot
#   fixes: the evaluator plots its own geometry, element-panel phases
#   wrap, and the polar dial draws the mirrored half of the pattern.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from array_common import ArrayCommon
from array_evaluator import ArrayEvaluator
from array_calculator import ArrayCalculator
from hamming_array import HammingArray
from uniform_array import UniformArray
from test_helpers import make_args


def _brute_directivity(amps, positions, phases):
    th = np.linspace(1e-6, np.pi - 1e-6, 200001)
    w = np.asarray(amps) * np.exp(1j * np.asarray(phases))
    af = np.abs((w[None, :] * np.exp(1j * 2 * np.pi * np.asarray(positions)[None, :]
                                     * np.cos(th)[:, None])).sum(1)) ** 2
    return 2 * af.max() / np.trapezoid(af * np.sin(th), th)


class TestDirectivity:
    @pytest.mark.parametrize("d", [0.5, 0.7, 1.2])
    def test_matches_brute_force_integral(self, d):
        h = HammingArray(make_args(elements=16, spacing=d))
        amps = h.amplitudes(16)
        pos = np.arange(16) * d
        assert h.directivity(16, d, amps) == pytest.approx(
            _brute_directivity(amps, pos, np.zeros(16)), rel=1e-3)

    def test_grating_lobes_no_longer_inflate_directivity(self):
        # Regression: eff * 2 N d / lambda reported 14.5 dB at d = 1.2 lambda;
        # the grating lobe drops the true value to ~9.7 dB.
        h = HammingArray(make_args(elements=16, spacing=1.2))
        d0_db = 10 * np.log10(h.directivity(16, 1.2, h.amplitudes(16)))
        assert d0_db == pytest.approx(9.73, abs=0.1)

    def test_uniform_broadside_half_wave_is_N(self):
        u = UniformArray(make_args(elements=10))
        assert u.directivity(10, 0.5) == pytest.approx(10.0, rel=1e-6)

    def test_steering_is_accounted_for(self):
        # Endfire (scan 0) directivity is about twice broadside at d = lambda/4.
        bs = UniformArray(make_args(elements=20, spacing=0.25, scan=90.0))
        ef = UniformArray(make_args(elements=20, spacing=0.25, scan=0.0))
        ratio = ef.directivity(20, 0.25) / bs.directivity(20, 0.25)
        assert 1.7 < ratio < 2.3

    def test_binomial_off_half_wave_uses_integral(self):
        shell = ArrayCalculator(["binomial_array", "-N", "6", "-d", "0.8",
                                 "--variable_return"])
        shell.main(shell.getArgs())
        _, _, d0_db = shell.getCalcedParams()
        from math import comb
        a = np.array([comb(5, n) for n in range(6)], float)
        brute = 10 * np.log10(_brute_directivity(a, np.arange(6) * 0.8, np.zeros(6)))
        assert d0_db == pytest.approx(brute, abs=0.02)


class TestHPBW:
    def test_endfire_full_width(self):
        # Regression: the walk stopped at theta = 0 and reported only half
        # of an endfire beam. Endfire is symmetric about the axis, so the
        # width is twice the one-sided half-power angle.
        u = UniformArray(make_args(elements=10, spacing=0.25, scan=0.0))
        amps = u.amplitudes(10)
        theta, af, _ = u.pattern_sweep(amps, n_points=40001)
        one_side = theta[np.argmax(af < 1 / np.sqrt(2))]
        assert u.hpbw_deg(amps) == pytest.approx(2 * one_side, abs=0.02)

    def test_broadside_unchanged(self):
        u = UniformArray(make_args(elements=10))
        # Balanis: uniform N=10, d=lambda/2 broadside HPBW ~ 10.2 deg
        assert u.hpbw_deg(u.amplitudes(10)) == pytest.approx(10.2, abs=0.1)

    def test_single_element_has_no_hpbw(self):
        c = ArrayCommon(make_args())
        assert c.hpbw_deg([1.0]) is None


class TestPeakSidelobe:
    def test_grating_lobe_counts(self):
        # Regression: lobes above -0.5 dB were ignored, so a uniform array at
        # d = 1.0 lambda (grating lobes at 0 and 180 deg) scored -13 dB.
        ev = ArrayEvaluator(make_args())
        res = ev.evaluate(np.arange(8) * 1.0)
        assert res["peak_sidelobe_db"] > -0.5

    def test_endpoint_lobe_counts(self):
        # A lobe sitting on theta = 0 must be seen even though the sweep
        # starts there (no interior local maximum).
        theta = np.linspace(0, 180, 181)
        af_db = np.full_like(theta, -40.0)
        af_db[80:101] = -np.abs(theta[80:101] - 90) / 2   # main beam at 90
        af_db[:5] = [-8, -9, -12, -20, -30]              # lobe peaking at 0
        assert ArrayCommon.peak_sidelobe(theta, af_db) == pytest.approx(-8.0)

    def test_normal_array_unchanged(self):
        ev = ArrayEvaluator(make_args())
        res = ev.evaluate(np.arange(16) * 0.5)
        assert res["peak_sidelobe_db"] == pytest.approx(-13.2, abs=0.2)


class TestPlots:
    def _capture(self, monkeypatch):
        seen = {}
        orig = ArrayCommon._rect_axes

        def spy(self, ax, theta, af_db, sll_db, floor):
            seen["theta"], seen["af_db"] = theta, af_db
            return orig(self, ax, theta, af_db, sll_db, floor)
        monkeypatch.setattr(ArrayCommon, "_rect_axes", spy)
        return seen

    def test_evaluator_plots_its_own_geometry(self, monkeypatch, tmp_path):
        # Regression: plot_pattern swept an equally spaced lambda/2 array
        # instead of the CSV geometry (up to ~52 dB off on the example file).
        pytest.importorskip("matplotlib")
        g = tmp_path / "g.csv"
        g.write_text("position_lambda\n0\n0.55\n1.15\n1.85\n2.65\n3.55\n4.55\n5.65\n")
        seen = self._capture(monkeypatch)
        shell = ArrayCalculator(["evaluate", "-g", str(g), "--save",
                                 str(tmp_path / "e.png"), "--plot-style", "rect"])
        shell.main(shell.getArgs())
        pos = np.loadtxt(g, skiprows=1)
        ev = ArrayEvaluator(make_args())
        _, _, truth = ev.pattern_sweep(np.ones(8), positions=pos, phases=np.zeros(8))
        np.testing.assert_allclose(seen["af_db"], np.clip(truth, -60, 0), atol=1e-9)

    def test_polar_draws_mirrored_half(self, tmp_path):
        mpl = pytest.importorskip("matplotlib")
        mpl.use("Agg")
        import matplotlib.pyplot as plt
        u = UniformArray(make_args(elements=8))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        theta, _, af_db = u.pattern_sweep(u.amplitudes(8))
        u._polar_axes(ax, theta, np.clip(af_db, -60, 0), None, -60)
        lo, hi = np.degrees(ax.get_xlim())
        assert hi - lo == pytest.approx(360.0)        # full dial, both halves
        xs = [ln.get_xdata() for ln in ax.get_lines()]
        assert any(np.min(x) < -1.0 for x in xs)      # mirrored trace present
        plt.close(fig)

    def test_element_phases_wrap(self, monkeypatch):
        mpl = pytest.importorskip("matplotlib")
        mpl.use("Agg")
        import matplotlib.pyplot as plt
        u = UniformArray(make_args(elements=8))
        fig, ax = plt.subplots()
        phases = np.radians([0, 170, 190, 360, 540, -200, 725, 90])
        u._element_axes(ax, np.arange(8) * 0.5, np.ones(8), phases)
        colors = ax.collections[-1].get_array()
        np.testing.assert_allclose(colors, [0, 170, -170, 0, -180, 160, 5, 90],
                                   atol=1e-9)
        plt.close(fig)


def test_array_factor_chunking_matches_direct_sum():
    # The chunked array factor must equal the one-shot matrix sum.
    rng = np.random.default_rng(0)
    N = 300
    amps = rng.uniform(0.2, 1.0, N)
    pos = np.sort(rng.uniform(0, 150, N))
    ph = rng.uniform(-np.pi, np.pi, N)
    theta = np.linspace(0, 180, 9001)
    c = ArrayCommon(make_args())
    got = c.array_factor(amps, theta, positions=pos, phases=ph)
    t = np.radians(theta)
    w = amps * np.exp(1j * ph)
    ref = np.abs((w[None, :] * np.exp(1j * 2 * np.pi * pos[None, :]
                                      * np.cos(t)[:, None])).sum(1))
    np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-9)
