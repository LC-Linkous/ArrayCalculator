#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_triangular_array.py
#
#   Tests for the triangular array (triangular_array.py): a linear taper
#   whose pattern is the square of the uniform array's, giving a first
#   sidelobe near -26 dB with decaying tails.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from triangular_array import TriangularArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


class TestTriangularSynthesis:
    def test_amplitudes_symmetric(self):
        t = TriangularArray(make_args(elements=10))
        for N in range(3, 16):
            a = t.amplitudes(N)
            assert np.allclose(a, a[::-1])

    def test_amplitudes_peak_at_center(self):
        t = TriangularArray(make_args(elements=11))
        a = t.amplitudes(11)
        assert np.argmax(a) == 5  # center element of an odd array

    def test_amplitudes_linear_ramp(self):
        # Even N: normalized amps should be an evenly spaced ramp up to 1.
        t = TriangularArray(make_args(elements=10))
        a = t.amplitudes(10)
        a = a / a.max()
        expected = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 0.8, 0.6, 0.4, 0.2])
        assert np.allclose(a, expected, atol=1e-9)


class TestTriangularPattern:
    @pytest.mark.parametrize("N", [10, 16, 20, 30])
    def test_sidelobe_near_26db(self, N):
        # First sidelobe of a triangular taper sits at about -26 dB.
        t = TriangularArray(make_args(elements=N))
        amps = t.amplitudes(N)
        assert peak_sidelobe_db(t, amps) == pytest.approx(-26.5, abs=1.0)

    def test_lower_sidelobes_than_uniform(self):
        from uniform_array import UniformArray
        N = 20
        t = TriangularArray(make_args(elements=N))
        u = UniformArray(make_args(elements=N))
        assert peak_sidelobe_db(t, t.amplitudes(N)) < peak_sidelobe_db(u, u.amplitudes(N))

    def test_main_beam_broadside(self):
        t = TriangularArray(make_args(elements=12))
        theta, af_norm, _ = t.pattern_sweep(t.amplitudes(12))
        assert theta[np.argmax(af_norm)] == pytest.approx(90.0, abs=0.5)


class TestTriangularCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["triangular_array", "-N", "12", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert len(amps) == 12
        assert amps.max() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
