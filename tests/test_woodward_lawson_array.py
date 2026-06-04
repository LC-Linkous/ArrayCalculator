#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_woodward_lawson_array.py
#
#   Tests for Woodward-Lawson shaped-beam synthesis
#   (woodward_lawson_array.py). Unlike the pencil-beam methods, the test
#   is that the REALIZED pattern approximates the requested SHAPE: a
#   flat-top holds a sector and falls off outside it; cosecant-squared
#   produces a broad, asymmetric coverage beam rather than a pencil.
##--------------------------------------------------------------------\

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


class TestCosecantSquared:
    def test_is_shaped_not_pencil(self):
        a = make_args(elements=24)
        a.shape, a.sector, a.floor, a.scan = "cosecant_squared", 30.0, 0.0, 90.0
        w = WoodwardLawsonArray(a)
        theta, af, _ = _realized(w, w.amplitudes(24))
        # a shaped beam has a wide -3 dB span, far wider than a pencil beam
        span = theta[af > 1 / np.sqrt(2)]
        assert (span.max() - span.min()) > 20.0


class TestWoodwardLawsonCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["woodward_lawson_array", "-N", "20",
                                 "--shape", "flat_top", "--variable_return"])
        shell.main(shell.getArgs())
        amps, phase_deg, d0_db = shell.getCalcedParams()
        assert len(amps) == 20


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
