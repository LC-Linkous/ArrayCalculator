#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_hann_array.py
#
#   Tests for the Hann array (hann_array.py): first sidelobe near
#   -31.4 dB, computed from the actual array factor.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from hann_array import HannArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


class TestHannArraySynthesis:
    def test_amplitudes_symmetric(self):
        obj = HannArray(make_args(elements=10))
        for N in range(3, 16):
            a = obj.amplitudes(N)
            assert np.allclose(a, a[::-1])

    def test_amplitudes_peak_at_center(self):
        obj = HannArray(make_args(elements=11))
        a = obj.amplitudes(11)
        assert np.argmax(a) == 5

    def test_amplitudes_positive(self):
        obj = HannArray(make_args(elements=12))
        for N in [5, 6, 8, 12, 20]:
            assert np.all(obj.amplitudes(N) >= 0)


class TestHannArrayPattern:
    @pytest.mark.parametrize("N", [10, 16, 20, 30])
    def test_sidelobe_signature(self, N):
        obj = HannArray(make_args(elements=N))
        amps = obj.amplitudes(N)
        assert peak_sidelobe_db(obj, amps) == pytest.approx(-31.4, abs=1)

    def test_main_beam_broadside(self):
        obj = HannArray(make_args(elements=12))
        theta, af_norm, _ = obj.pattern_sweep(obj.amplitudes(12))
        assert theta[np.argmax(af_norm)] == pytest.approx(90.0, abs=0.5)


class TestHannArrayCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["hann_array", "-N", "12", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert len(amps) == 12
        assert amps.max() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


class TestHannCosineSquaredEquivalence:
    def test_identical_to_cosine_squared(self):
        # Hann and cosine-squared are the same distribution
        # (0.5 + 0.5 cos(2x) == cos^2(x)); verify byte-for-byte across N.
        from cosine_squared_array import CosineSquaredArray
        h = HannArray(make_args(elements=10))
        c = CosineSquaredArray(make_args(elements=10))
        for N in [5, 6, 8, 11, 16, 20]:
            assert np.allclose(h.amplitudes(N), c.amplitudes(N), atol=1e-12)
