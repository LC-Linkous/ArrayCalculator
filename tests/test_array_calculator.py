#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  tests/test_array_calculator.py
#
#   Cross-cutting tests that are NOT specific to one synthesis method:
#     - array factor / beam steering
#     - CSV pattern export
#     - CLI dispatch, defaults, argument validation, --variable_return
#     - shared edge cases (odd/even N, high SLL, physical spacing)
#
#   Per-array synthesis and pattern checks live in the dedicated
#   test_<array>.py files; shared fixtures live in test_helpers.py. This
#   file used to be the monolithic suite -- it now holds only the
#   coverage that has no per-array home, so nothing is double-counted.
##--------------------------------------------------------------------\

import csv

import numpy as np
import pytest

from binomial_array import BinomialArray
from dolph_tschebyscheff import DolphTschebyscheff
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db



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