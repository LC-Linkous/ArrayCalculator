#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_blackman_array.py
#
#   Tests for the Blackman array (blackman_array.py): first sidelobe near
#   -57.5 dB, computed from the actual array factor.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from blackman_array import BlackmanArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


class TestBlackmanArraySynthesis:
    def test_amplitudes_symmetric(self):
        obj = BlackmanArray(make_args(elements=10))
        for N in range(3, 16):
            a = obj.amplitudes(N)
            assert np.allclose(a, a[::-1])

    def test_amplitudes_peak_at_center(self):
        obj = BlackmanArray(make_args(elements=11))
        a = obj.amplitudes(11)
        assert np.argmax(a) == 5

    def test_amplitudes_positive(self):
        obj = BlackmanArray(make_args(elements=12))
        for N in [5, 6, 8, 12, 20]:
            assert np.all(obj.amplitudes(N) >= 0)


class TestBlackmanArrayPattern:
    @pytest.mark.parametrize("N", [10, 16, 20, 30])
    def test_sidelobe_signature(self, N):
        obj = BlackmanArray(make_args(elements=N))
        amps = obj.amplitudes(N)
        assert peak_sidelobe_db(obj, amps) == pytest.approx(-57.5, abs=2.5)

    def test_main_beam_broadside(self):
        obj = BlackmanArray(make_args(elements=12))
        theta, af_norm, _ = obj.pattern_sweep(obj.amplitudes(12))
        assert theta[np.argmax(af_norm)] == pytest.approx(90.0, abs=0.5)


class TestBlackmanArrayCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["blackman_array", "-N", "12", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert len(amps) == 12
        assert amps.max() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
