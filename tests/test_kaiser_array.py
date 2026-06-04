#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_kaiser_array.py
#
#   Tests for the Kaiser array (kaiser_array.py): parametric taper whose
#   sidelobe level falls monotonically with the shape parameter beta.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from kaiser_array import KaiserArray, _i0
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


class TestKaiserBessel:
    @pytest.mark.parametrize("x,ref", [
        (0.0, 1.0), (1.0, 1.2660658), (2.0, 2.2795853), (5.0, 27.2398718),
    ])
    def test_i0_matches_reference(self, x, ref):
        assert float(_i0(np.array([x]))[0]) == pytest.approx(ref, rel=1e-6)


class TestKaiserSynthesis:
    def test_amplitudes_symmetric(self):
        k = KaiserArray(make_args(elements=10))
        for N in range(3, 16):
            a = k.amplitudes(N, beta=6.0)
            assert np.allclose(a, a[::-1])

    def test_peak_at_center(self):
        k = KaiserArray(make_args(elements=11))
        assert np.argmax(k.amplitudes(11, beta=6.0)) == 5

    def test_beta_zero_is_uniform(self):
        # beta = 0 -> I0(0)=1 everywhere -> uniform excitation.
        k = KaiserArray(make_args(elements=12))
        assert np.allclose(k.amplitudes(12, beta=0.0), 1.0)


class TestKaiserPattern:
    def test_sidelobes_drop_with_beta(self):
        k = KaiserArray(make_args(elements=32))
        levels = [peak_sidelobe_db(k, k.amplitudes(32, beta=b)) for b in [2, 4, 6, 8]]
        # strictly decreasing sidelobe level as beta grows
        assert all(levels[i] > levels[i + 1] for i in range(len(levels) - 1))

    def test_main_beam_broadside(self):
        k = KaiserArray(make_args(elements=16))
        theta, af, _ = k.pattern_sweep(k.amplitudes(16, beta=6.0))
        assert theta[np.argmax(af)] == pytest.approx(90.0, abs=0.5)


class TestKaiserCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["kaiser_array", "-N", "16", "-beta", "8",
                                 "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert len(amps) == 16
        assert amps.max() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
