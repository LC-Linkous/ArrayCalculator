#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_villeneuve_array.py
#
#   Tests for the Villeneuve array (villeneuve_array.py): discrete analog
#   of Taylor. Defining guarantee: the peak sidelobe lands at the design
#   level, with near-in sidelobes held there and the rest decaying.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from villeneuve_array import VilleneuveArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


class TestVilleneuveSynthesis:
    def test_length_and_symmetry(self):
        v = VilleneuveArray(make_args(elements=20, sidelobe_level=30, nbar=6))
        for N in [16, 20, 24, 30]:
            amps, _, _ = v.amplitudes(N, 30, 6)
            assert len(amps) == N
            assert np.allclose(amps, amps[::-1], atol=1e-6)

    def test_sigma_geq_one(self):
        v = VilleneuveArray(make_args(elements=20, sidelobe_level=30, nbar=6))
        _, _, sigma = v.amplitudes(20, 30, 6)
        assert sigma >= 1.0


class TestVilleneuvePattern:
    @pytest.mark.parametrize("N,sll,nbar", [
        (20, 30, 6), (20, 25, 5), (30, 35, 8), (40, 40, 10),
    ])
    def test_peak_sidelobe_hits_target(self, N, sll, nbar):
        v = VilleneuveArray(make_args(elements=N, sidelobe_level=sll, nbar=nbar))
        amps, _, _ = v.amplitudes(N, sll, nbar)
        assert peak_sidelobe_db(v, amps) == pytest.approx(-sll, abs=0.5)

    def test_main_beam_broadside(self):
        v = VilleneuveArray(make_args(elements=20, sidelobe_level=30, nbar=6))
        amps, _, _ = v.amplitudes(20, 30, 6)
        theta, af, _ = v.pattern_sweep(amps)
        assert theta[np.argmax(af)] == pytest.approx(90.0, abs=0.5)


class TestVilleneuveCLI:
    def test_variable_return(self):
        shell = ArrayCalculator(["villeneuve_array", "-N", "20", "-sll", "30",
                                 "-nbar", "6", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert len(amps) == 20
        assert amps.max() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
