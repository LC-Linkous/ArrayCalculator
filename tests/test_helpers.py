#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_array_calculator.py
#
#   Test suite for the array calculator. Run with:  pytest -v
#
#   Coverage:
#     - Binomial synthesis vs. the EGRE 540 notes (Examples 1-3, Problems)
#     - Dolph-Tschebyscheff synthesis vs. notes (Examples 4-6)
#     - Physical ground truth: computed patterns hit the target sidelobe
#       level, binomial arrays have no sidelobes, etc.
#     - Array-factor / pattern-sweep behavior (broadside + steered)
#     - CSV export and --variable_return
#     - CLI dispatch and argument handling
##--------------------------------------------------------------------\

import csv
import math
from types import SimpleNamespace

import numpy as np
import pytest

from binomial_array import BinomialArray
from dolph_tschebyscheff import DolphTschebyscheff
from array_calculator import ArrayCalculator


# --- helpers -----------------------------------------------------------
def make_args(**kwargs):
    """Build an args namespace with the defaults the classes expect."""
    defaults = dict(
        elements=None, frequency=None, spacing=0.5, scan=90.0,
        unit="centimeter", csv=None, plot=False, verbose=False,
        variable_return=False, sidelobe_level=None, norm="edge",
    )
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


def peak_sidelobe_db(obj, amps):
    """Highest sidelobe level (dB) of the array's actual pattern."""
    theta, af_norm, af_db = obj.pattern_sweep(amps, n_points=20001)
    peaks = [af_db[i] for i in range(1, len(af_db) - 1)
             if af_db[i] >= af_db[i - 1] and af_db[i] >= af_db[i + 1]
             and af_db[i] < -0.5]
    return max(peaks) if peaks else -np.inf


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


# ======================================================================
# ARRAY FACTOR / STEERING
# ======================================================================
class TestArrayFactor:
    def test_broadside_peak_at_90(self):
        b = BinomialArray(make_args(elements=6, scan=90.0))
        amps = b.amplitudes(6)
        theta, af, _ = b.pattern_sweep(amps)
        assert theta[np.argmax(af)] == pytest.approx(90.0, abs=0.5)

    def test_steered_peak_moves(self):
        # Steering to 60 deg should move the main beam toward 60 deg.
        dt = DolphTschebyscheff(make_args(elements=8, sidelobe_level=25, scan=60.0))
        amps, _, _ = dt.amplitudes(8, 25)
        theta, af, _ = dt.pattern_sweep(amps)
        peak_angle = theta[np.argmax(af)]
        assert peak_angle == pytest.approx(60.0, abs=2.0)

    def test_uniform_array_factor_known_value(self):
        # Uniform 2-element, d=lambda/2, broadside: AF at endfire (theta=0)
        # = |1 + exp(j*pi*cos0)| with broadside phasing = |1 + e^{j*pi}| = 0.
        b = BinomialArray(make_args(elements=2))
        af = b.array_factor([1.0, 1.0], [0.0], d_over_lambda=0.5, scan_deg=90.0)
        assert af[0] == pytest.approx(0.0, abs=1e-9)
        af90 = b.array_factor([1.0, 1.0], [90.0], d_over_lambda=0.5, scan_deg=90.0)
        assert af90[0] == pytest.approx(2.0, abs=1e-9)

    def test_pattern_normalized_to_one(self):
        dt = DolphTschebyscheff(make_args(elements=10, sidelobe_level=26))
        amps, _, _ = dt.amplitudes(10, 26)
        _, af_norm, af_db = dt.pattern_sweep(amps)
        assert af_norm.max() == pytest.approx(1.0)
        assert af_db.max() == pytest.approx(0.0, abs=1e-6)


# ======================================================================
# CSV EXPORT
# ======================================================================
class TestCsvExport:
    def test_csv_written_with_headers(self, tmp_path):
        out = tmp_path / "pattern.csv"
        b = BinomialArray(make_args(elements=6))
        amps = b.amplitudes(6)
        b.export_pattern_csv(str(out), amps)
        assert out.exists()
        with open(out) as f:
            rows = list(csv.reader(f))
        assert rows[0] == ["theta_deg", "AF_linear", "AF_dB"]
        assert len(rows) > 100  # header + sweep points

    def test_csv_values_parseable(self, tmp_path):
        out = tmp_path / "p.csv"
        dt = DolphTschebyscheff(make_args(elements=8, sidelobe_level=25))
        amps, _, _ = dt.amplitudes(8, 25)
        dt.export_pattern_csv(str(out), amps)
        with open(out) as f:
            rows = list(csv.DictReader(f))
        thetas = [float(r["theta_deg"]) for r in rows]
        lin = [float(r["AF_linear"]) for r in rows]
        assert thetas[0] == pytest.approx(0.0)
        assert thetas[-1] == pytest.approx(180.0)
        assert max(lin) == pytest.approx(1.0, abs=1e-6)


# ======================================================================
# CLI / DISPATCH
# ======================================================================
class TestCLI:
    def test_binomial_variable_return(self):
        shell = ArrayCalculator(["binomial_array", "-N", "6", "--variable_return"])
        shell.main(shell.getArgs())
        amps, hpbw, d0_db = shell.getCalcedParams()
        assert np.allclose(amps, [0.1, 0.5, 1.0, 1.0, 0.5, 0.1])
        assert hpbw == pytest.approx(27.2, abs=0.1)
        assert d0_db == pytest.approx(6.37, abs=0.02)

    def test_dolph_variable_return(self):
        shell = ArrayCalculator(["dolph_tschebyscheff", "-N", "9", "-sll", "25",
                                 "--norm", "center", "--variable_return"])
        shell.main(shell.getArgs())
        amps, R, z0, d0_db = shell.getCalcedParams()
        assert amps.max() == pytest.approx(1.0)        # center-normalized
        assert R == pytest.approx(17.78, abs=0.01)
        assert z0 == pytest.approx(1.101, abs=0.001)

    def test_dispatch_sets_subparser_name(self):
        shell = ArrayCalculator(["binomial_array", "-N", "4"])
        assert shell.getArgs().subparser_name == "binomial_array"

    def test_defaults(self):
        shell = ArrayCalculator(["binomial_array", "-N", "5"])
        args = shell.getArgs()
        assert args.spacing == 0.5
        assert args.scan == 90.0
        assert args.variable_return is False

    def test_frequency_optional(self):
        # No -f provided: should still run and return values.
        shell = ArrayCalculator(["binomial_array", "-N", "5", "--variable_return"])
        shell.main(shell.getArgs())
        assert shell.getCalcedParams() is not None

    def test_missing_required_elements_errors(self):
        with pytest.raises(SystemExit):
            ArrayCalculator(["binomial_array"])

    def test_dolph_missing_sll_errors(self):
        with pytest.raises(SystemExit):
            ArrayCalculator(["dolph_tschebyscheff", "-N", "10"])

    def test_invalid_norm_choice_errors(self):
        with pytest.raises(SystemExit):
            ArrayCalculator(["dolph_tschebyscheff", "-N", "10", "-sll", "26",
                             "--norm", "middle"])


# ======================================================================
# EDGE CASES
# ======================================================================
class TestEdgeCases:
    def test_odd_and_even_N_both_work(self):
        dt = DolphTschebyscheff(make_args(elements=7, sidelobe_level=25))
        for N in [5, 6, 7, 8, 9, 10, 11, 12]:
            amps, _, _ = dt.amplitudes(N, 25)
            assert len(amps) == N
            assert np.all(amps > 0)

    def test_high_sidelobe_level(self):
        dt = DolphTschebyscheff(make_args(elements=10, sidelobe_level=40))
        amps, _, _ = dt.amplitudes(10, 40)
        assert peak_sidelobe_db(dt, amps) == pytest.approx(-40, abs=0.2)

    def test_spacing_fraction_changes_physical_d(self):
        # With frequency set, d should scale with the spacing fraction.
        lam = 3e8 / 3e9  # 0.1 m at 3 GHz
        b = BinomialArray(make_args(elements=6, frequency=3e9, spacing=0.4))
        assert b.wavelength() == pytest.approx(lam)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))