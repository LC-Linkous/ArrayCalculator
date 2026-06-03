#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_dolph_tschebyscheff.py
#
#   Tests for the Dolph-Tschebyscheff synthesis (dolph_tschebyscheff.py),
#   checked against the EGRE 540 notes (Examples 4-6) and the defining
#   physical guarantee: equiripple sidelobes at exactly the specified
#   level, with a narrower beam than binomial for equal N.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from binomial_array import BinomialArray
from dolph_tschebyscheff import DolphTschebyscheff
from test_helpers import make_args, peak_sidelobe_db


# ======================================================================
# DOLPH-TSCHEBYSCHEFF SYNTHESIS
# ======================================================================
class TestDolphSynthesis:
    def test_z0_example4(self):
        # Notes Example 4: N=10, 26 dB -> z0 = 1.085
        dt = DolphTschebyscheff(make_args(elements=10, sidelobe_level=26))
        R = dt.voltage_ratio(26)
        assert R == pytest.approx(19.953, abs=0.01)
        assert dt.z0(10, R) == pytest.approx(1.085, abs=0.001)

    def test_z0_example5(self):
        # Notes Example 5: N=6, 30 dB -> z0 = 1.358
        dt = DolphTschebyscheff(make_args(elements=6, sidelobe_level=30))
        R = dt.voltage_ratio(30)
        assert R == pytest.approx(31.62, abs=0.01)
        assert dt.z0(6, R) == pytest.approx(1.358, abs=0.01)

    def test_z0_example6(self):
        # Notes Example 6: N=9, 25 dB -> z0 = 1.101
        dt = DolphTschebyscheff(make_args(elements=9, sidelobe_level=25))
        R = dt.voltage_ratio(25)
        assert R == pytest.approx(17.78, abs=0.01)
        assert dt.z0(9, R) == pytest.approx(1.101, abs=0.001)

    def test_voltage_ratio_known_points(self):
        dt = DolphTschebyscheff(make_args(elements=6, sidelobe_level=20))
        assert dt.voltage_ratio(20) == pytest.approx(10.0, abs=1e-9)
        assert dt.voltage_ratio(0) == pytest.approx(1.0, abs=1e-9)

    def test_coeffs_symmetric(self):
        dt = DolphTschebyscheff(make_args(elements=10, sidelobe_level=26))
        for N in [6, 7, 8, 9, 10, 11]:
            amps, _, _ = dt.amplitudes(N, 25)
            assert np.allclose(amps, amps[::-1], atol=1e-6)

    def test_edge_vs_center_normalization(self):
        # Edge-norm divides by the smallest element; center-norm by the largest.
        dt = DolphTschebyscheff(make_args(elements=10, sidelobe_level=26))
        amps, _, _ = dt.amplitudes(10, 26)
        edge = amps / amps.min()
        center = amps / amps.max()
        assert edge.min() == pytest.approx(1.0)
        assert center.max() == pytest.approx(1.0)
        # the two are scalar multiples of each other
        ratio = edge / center
        assert np.allclose(ratio, ratio[0])

    def test_coeffs_match_corrected_example4(self):
        # Corrected Example 4 values (verified to hit -26 dB), edge-normalized.
        dt = DolphTschebyscheff(make_args(elements=10, sidelobe_level=26))
        amps, _, _ = dt.amplitudes(10, 26)
        edge = amps / amps.min()
        expected = [1.0, 1.355, 1.968, 2.479, 2.769]
        assert np.allclose(edge[:5], expected, atol=0.01)


# ======================================================================
# DOLPH PHYSICAL PATTERN  (the real ground truth)
# ======================================================================
class TestDolphPattern:
    @pytest.mark.parametrize("N,sll", [
        (10, 26), (6, 30), (9, 25), (8, 25), (12, 20), (7, 30), (5, 20),
    ])
    def test_sidelobe_level_hits_target(self, N, sll):
        # The defining guarantee: actual peak sidelobe == specified SLL.
        dt = DolphTschebyscheff(make_args(elements=N, sidelobe_level=sll))
        amps, _, _ = dt.amplitudes(N, sll)
        assert peak_sidelobe_db(dt, amps) == pytest.approx(-sll, abs=0.1)

    def test_equiripple_sidelobes(self):
        # All sidelobes should be at essentially the same level (equiripple).
        dt = DolphTschebyscheff(make_args(elements=12, sidelobe_level=25))
        amps, _, _ = dt.amplitudes(12, 25)
        _, _, af_db = dt.pattern_sweep(amps, n_points=20001)
        peaks = [af_db[i] for i in range(1, len(af_db) - 1)
                 if af_db[i] >= af_db[i - 1] and af_db[i] >= af_db[i + 1]
                 and af_db[i] < -0.5]
        assert max(peaks) - min(peaks) < 0.2  # flat to within 0.2 dB

    def test_narrower_beam_than_binomial(self):
        # For equal N, Dolph should have a narrower main beam than binomial.
        N = 10
        b = BinomialArray(make_args(elements=N))
        dt = DolphTschebyscheff(make_args(elements=N, sidelobe_level=26))
        bamp = b.amplitudes(N)
        damp, _, _ = dt.amplitudes(N, 26)

        def first_null_width(obj, amps):
            theta, af, _ = obj.pattern_sweep(amps, n_points=40001)
            peak = np.argmax(af)
            # walk right from the peak to first local min
            i = peak
            while i + 1 < len(af) and af[i + 1] <= af[i]:
                i += 1
            return theta[i] - theta[peak]

        assert first_null_width(dt, damp) < first_null_width(b, bamp)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))