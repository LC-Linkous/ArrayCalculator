#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_uniform_array.py
#
#   Tests for the uniform array (uniform_array.py), the reference case:
#   equal excitation, directivity D0 = N at d = lambda/2, and the classic
#   ~-13.2 dB first sidelobe.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from uniform_array import UniformArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


# ======================================================================
# UNIFORM ARRAY
# ======================================================================
class TestUniformArray:
    def test_amplitudes_all_equal(self):
        u = UniformArray(make_args(elements=8))
        assert u.amplitudes(8).tolist() == [1.0] * 8

    def test_directivity_equals_N_at_half_wave(self):
        # D0 = N for a uniform broadside array at d = lambda/2.
        u = UniformArray(make_args(elements=10))
        assert u.directivity(10, 0.5) == pytest.approx(10.0)
        assert u.directivity(20, 0.5) == pytest.approx(20.0)

    def test_sidelobe_level_about_13db(self):
        # Uniform array's first sidelobe is the classic ~-13.2 dB.
        u = UniformArray(make_args(elements=20))
        amps = u.amplitudes(20)
        assert peak_sidelobe_db(u, amps) == pytest.approx(-13.2, abs=0.3)

    def test_variable_return(self):
        shell = ArrayCalculator(["uniform_array", "-N", "10", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert np.allclose(amps, 1.0)
        assert d0_db == pytest.approx(10.0, abs=0.01)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))