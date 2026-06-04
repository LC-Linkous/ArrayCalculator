#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_bartlett_array.py
#
#   Tests for the Bartlett array (bartlett_array.py): zero-endpoint
#   triangular taper, first sidelobe near -26 dB.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from bartlett_array import BartlettArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


class TestBartlettSynthesis:
    def test_amplitudes_symmetric(self):
        b = BartlettArray(make_args(elements=10))
        for N in range(3, 16):
            a = b.amplitudes(N)
            assert np.allclose(a, a[::-1])

    def test_edges_are_zero(self):
        # The defining feature vs triangular_array: endpoints go to zero.
        b = BartlettArray(make_args(elements=11))
        a = b.amplitudes(11)
        assert a[0] == pytest.approx(0.0)
        assert a[-1] == pytest.approx(0.0)

    def test_peak_at_center(self):
        b = BartlettArray(make_args(elements=11))
        assert np.argmax(b.amplitudes(11)) == 5


class TestBartlettPattern:
    @pytest.mark.parametrize("N", [16, 20, 32])
    def test_sidelobe_near_26db(self, N):
        b = BartlettArray(make_args(elements=N))
        assert peak_sidelobe_db(b, b.amplitudes(N)) == pytest.approx(-26.4, abs=1.0)

    def test_main_beam_broadside(self):
        b = BartlettArray(make_args(elements=12))
        theta, af, _ = b.pattern_sweep(b.amplitudes(12))
        assert theta[np.argmax(af)] == pytest.approx(90.0, abs=0.5)


class TestBartlettCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["bartlett_array", "-N", "12", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert len(amps) == 12
        assert amps.max() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
