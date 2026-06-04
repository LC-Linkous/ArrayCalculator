#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_binomial_array.py
#
#   Tests for the binomial array synthesis (binomial_array.py), checked
#   against the EGRE 540 notes (Examples 1-3, Problem 2) and the
#   defining physical property: no sidelobes.
##--------------------------------------------------------------------\

import math

import numpy as np
import pytest

from binomial_array import BinomialArray
from test_helpers import make_args, peak_sidelobe_db


# ======================================================================
# BINOMIAL SYNTHESIS
# ======================================================================
class TestBinomialSynthesis:
    def test_amplitudes_pascals_triangle(self):
        b = BinomialArray(make_args(elements=6))
        assert b.amplitudes(2).tolist() == [1, 1]
        assert b.amplitudes(5).tolist() == [1, 4, 6, 4, 1]
        assert b.amplitudes(6).tolist() == [1, 5, 10, 10, 5, 1]
        assert b.amplitudes(9).tolist() == [1, 8, 28, 56, 70, 56, 28, 8, 1]
        assert b.amplitudes(10).tolist() == [1, 9, 36, 84, 126, 126, 84, 36, 9, 1]

    def test_amplitudes_symmetric(self):
        b = BinomialArray(make_args(elements=7))
        for N in range(2, 15):
            a = b.amplitudes(N)
            assert np.allclose(a, a[::-1])

    def test_example1_hpbw_directivity(self):
        # Notes Example 1: N=6 at 3 GHz -> HPBW 27.2 deg, D0 6.37 dB
        b = BinomialArray(make_args(elements=6))
        assert math.degrees(b.hpbw(6)) == pytest.approx(27.2, abs=0.1)
        D0 = b.directivity(6, 0.5)
        assert 10 * math.log10(D0) == pytest.approx(6.37, abs=0.02)

    def test_example2_hpbw(self):
        # Notes Example 2: N=9 -> HPBW 21.5 deg
        b = BinomialArray(make_args(elements=9))
        assert math.degrees(b.hpbw(9)) == pytest.approx(21.5, abs=0.1)

    def test_example3_hpbw_directivity(self):
        # Notes Example 3: N=10 -> HPBW 20.2 deg, D0 7.48 dB
        b = BinomialArray(make_args(elements=10))
        assert math.degrees(b.hpbw(10)) == pytest.approx(20.2, abs=0.1)
        assert 10 * math.log10(b.directivity(10, 0.5)) == pytest.approx(7.48, abs=0.02)

    def test_problem2_hpbw(self):
        # Notes Problem 2: N=10, d=lambda/2 -> HPBW 20.2 deg
        b = BinomialArray(make_args(elements=10))
        assert math.degrees(b.hpbw(10)) == pytest.approx(20.2, abs=0.1)

    def test_directivity_two_forms_agree(self):
        # D0 = 1.77*sqrt(1+2L/lambda) should equal 1.77*sqrt(N) at d=lambda/2
        b = BinomialArray(make_args(elements=6))
        for N in [5, 6, 7, 9, 10]:
            assert b.directivity(N, 0.5) == pytest.approx(1.77 * math.sqrt(N), rel=1e-9)

    def test_directivity_is_approximation_above_integral(self):
        # The closed form is a large-array fit and reads a few percent high
        # versus the directivity integrated from the pattern. Documented here
        # so the gap is a known property, not a surprise: it should be close
        # (within ~10%) and on the high side.
        b = BinomialArray(make_args(elements=6))
        theta = np.linspace(1e-6, np.pi - 1e-6, 400001)
        for N in [6, 9, 10]:
            amps = b.amplitudes(N)
            n = np.arange(N)
            psi = 2.0 * np.pi * 0.5 * np.cos(theta)[:, None]
            af = np.abs((amps[None, :] * np.exp(1j * n[None, :] * psi)).sum(axis=1))
            power = af ** 2
            integral = (2.0 * power[np.argmin(np.abs(theta - np.pi / 2))]
                        / np.trapezoid(power * np.sin(theta), theta))
            formula = b.directivity(N, 0.5)
            assert formula >= integral
            assert formula == pytest.approx(integral, rel=0.10)


# ======================================================================
# BINOMIAL PHYSICAL PATTERN
# ======================================================================
class TestBinomialPattern:
    def test_no_sidelobes(self):
        # The defining property of a binomial array: no sidelobes.
        b = BinomialArray(make_args(elements=8))
        amps = b.amplitudes(8) / b.amplitudes(8).max()
        assert peak_sidelobe_db(b, amps) < -40  # effectively none

    def test_main_beam_broadside(self):
        b = BinomialArray(make_args(elements=6))
        amps = b.amplitudes(6)
        theta, af_norm, _ = b.pattern_sweep(amps)
        # peak should be at theta = 90 deg (broadside)
        assert theta[np.argmax(af_norm)] == pytest.approx(90.0, abs=0.5)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))