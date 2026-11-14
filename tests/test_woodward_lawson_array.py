#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_woodward_lawson_array.py
#
#   Tests for Woodward-Lawson shaped-beam synthesis
#   (woodward_lawson_array.py). Unlike the pencil-beam methods, the test
#   is that the REALIZED pattern approximates the requested SHAPE: a
#   flat-top holds a sector and falls off outside it; cosecant-squared
#   produces a broad, one-sided coverage beam rather than a pencil. The
#   CLI tests check what is actually exported, not just the helper math.
##--------------------------------------------------------------------\

import csv

import numpy as np
import pytest

from woodward_lawson_array import WoodwardLawsonArray
from array_calculator import ArrayCalculator
from test_helpers import make_args


def _realized(obj, amps_c, n_points=2001):
    mag = np.abs(amps_c)
    ph = np.angle(amps_c)
    theta = np.linspace(0, 180, n_points)
    af = obj.array_factor(mag, theta, phases=ph)
    af = af / af.max()
    return theta, af, 20 * np.log10(np.clip(af, 1e-12, None))


class TestFlatTop:
    def test_holds_sector_and_drops_outside(self):
        a = make_args(elements=24)
        a.shape, a.sector, a.floor = "flat_top", 30.0, 0.0
        w = WoodwardLawsonArray(a)
        theta, af, afdb = _realized(w, w.amplitudes(24))
        # in-sector (60..120) should be high; far out (near 30) should be low
        insec = (theta >= 65) & (theta <= 115)
        farout = np.argmin(np.abs(theta - 35))
        assert afdb[insec].min() > -6.0          # roughly held up across the sector
        assert afdb[farout] < -15.0              # suppressed outside

    def test_broadside_excitations_real(self):
        a = make_args(elements=20)
        a.shape, a.sector, a.floor = "flat_top", 30.0, 0.0
        w = WoodwardLawsonArray(a)
        amps = w.amplitudes(20)
        assert np.max(np.abs(amps.imag)) < 1e-9  # symmetric target -> real


class TestSteering:
    def test_steered_sector_lands_on_scan_side(self):
        # Regression: the composing beams used exp(+j...), which put a
        # steered sector at 180 - theta (106..130 deg instead of 45..75).
        a = make_args(elements=24, scan=60.0)
        a.shape, a.sector, a.floor = "flat_top", 15.0, 0.0
        w = WoodwardLawsonArray(a)
        theta, af, afdb = _realized(w, w.amplitudes(24))
        band = theta[afdb > -3.0]
        assert 40.0 < band.min() < 55.0 and 70.0 < band.max() < 80.0


class TestCosecantSquared:
    def _w(self, scan=80.0, sector=40.0):
        a = make_args(elements=24, scan=scan)
        a.shape, a.sector, a.floor = "cosecant_squared", sector, 0.0
        return WoodwardLawsonArray(a)

    def test_is_shaped_not_pencil(self):
        # The realized pattern tracks the csc target across the fill region
        # (within the usual Woodward-Lawson ripple) instead of a pencil beam.
        w = self._w()
        theta, _, afdb = _realized(w, w.amplitudes(24))
        fill = (theta >= 45) & (theta <= 75)
        target_db = 20 * np.log10(w._target(theta[fill]))
        assert np.max(np.abs(afdb[fill] - target_db)) < 4.5
        span = theta[afdb > -15.0]
        assert (span.max() - span.min()) > 30.0

    def test_is_one_sided(self):
        # Regression: the old target used sin(theta), which is symmetric
        # about broadside. Coverage at 60 deg (inside the fill, above the
        # horizon) must far exceed its mirror image at 120 deg.
        w = self._w()
        theta, _, afdb = _realized(w, w.amplitudes(24))
        i = lambda t: np.argmin(np.abs(theta - t))
        assert afdb[i(60)] > afdb[i(120)] + 15.0

    def test_target_follows_csc_law(self):
        w = self._w(scan=80.0, sector=40.0)
        th = np.array([80.0, 60.0, 45.0, 100.0])
        t = w._target(th)
        assert t[0] == pytest.approx(1.0)
        assert t[1] == pytest.approx(np.cos(np.radians(80)) / np.cos(np.radians(60)))
        assert t[3] == 0.0                       # wrong side of broadside

    def test_looking_down_mirrors(self):
        up, down = self._w(scan=80.0), self._w(scan=100.0)
        th = np.linspace(0, 180, 181)
        np.testing.assert_allclose(up._target(th), down._target(180 - th))

    def test_empty_sector_is_a_clear_error(self):
        a = make_args(elements=4, scan=40.0)
        a.shape, a.sector, a.floor = "flat_top", 1.0, 0.0
        with pytest.raises(ValueError, match="sector"):
            WoodwardLawsonArray(a).amplitudes(4)

    def test_broadside_scan_rejected(self):
        with pytest.raises(ValueError):
            self._w(scan=90.0).amplitudes(24)


class TestWoodwardLawsonCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["woodward_lawson_array", "-N", "20",
                                 "--shape", "flat_top", "--variable_return"])
        shell.main(shell.getArgs())
        amps, phase_deg, d0_db = shell.getCalcedParams()
        assert len(amps) == 20

    def test_exported_pattern_is_the_shaped_beam(self, tmp_path):
        # Regression: the CSV used |a_n| only, dropping the 180 deg phases,
        # which turned the flat-top into a pencil beam (-11 dB at 70 deg).
        out = tmp_path / "wl.csv"
        shell = ArrayCalculator(["woodward_lawson_array", "-N", "24",
                                 "--shape", "flat_top", "--sector", "30",
                                 "--csv", str(out)])
        shell.main(shell.getArgs())
        rows = list(csv.DictReader(open(out)))
        theta = np.array([float(r["theta_deg"]) for r in rows])
        db = np.array([float(r["AF_dB"]) for r in rows])
        insec = (theta >= 65) & (theta <= 115)
        assert db[insec].min() > -6.0
        assert db[np.argmin(np.abs(theta - 35))] < -15.0

    def test_directivity_matches_realized_pattern(self):
        a = make_args(elements=24)
        a.shape, a.sector, a.floor = "flat_top", 30.0, 0.0
        w = WoodwardLawsonArray(a)
        amps = w.amplitudes(24)
        th = np.linspace(1e-6, np.pi - 1e-6, 200001)
        pos = np.arange(24) * 0.5
        af = np.abs((amps[None, :] * np.exp(1j * 2 * np.pi * pos[None, :]
                                            * np.cos(th)[:, None])).sum(1)) ** 2
        D = 2 * af.max() / np.trapezoid(af * np.sin(th), th)
        assert w.directivity(24, 0.5, amps) == pytest.approx(D, rel=1e-3)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
