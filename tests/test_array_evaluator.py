#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_array_evaluator.py
#
#   Tests for the evaluator (array_evaluator.py): the optimizer hook that
#   scores an arbitrary geometry. The defining guarantee is the
#   round-trip -- a uniform geometry scored by the evaluator must match
#   what the dedicated UniformArray produces -- plus correct response to
#   per-element phase (steering) and non-uniform spacing.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from array_evaluator import ArrayEvaluator
from uniform_array import UniformArray
from array_calculator import ArrayCalculator
from test_helpers import make_args


class TestEvaluatorRoundTrip:
    def test_uniform_geometry_matches_uniform_array(self):
        N, d = 10, 0.5
        ev = ArrayEvaluator(make_args())
        res = ev.evaluate(np.arange(N) * d, np.ones(N), np.zeros(N))
        u = UniformArray(make_args(elements=N))
        assert res["directivity"] == pytest.approx(u.directivity(N, d), rel=1e-3)
        assert res["peak_theta_deg"] == pytest.approx(90.0, abs=0.5)
        assert res["peak_sidelobe_db"] == pytest.approx(-13.0, abs=0.5)

    def test_defaults_uniform_amps_zero_phase(self):
        ev = ArrayEvaluator(make_args())
        res = ev.evaluate(np.arange(8) * 0.5)  # amps/phases default
        assert np.allclose(res["amplitudes"], 1.0)
        assert np.allclose(res["phases"], 0.0)


class TestEvaluatorPhaseSteering:
    def test_phase_steers_beam(self):
        N, d = 10, 0.5
        ev = ArrayEvaluator(make_args())
        pos = np.arange(N) * d
        beta = -2 * np.pi * pos * np.cos(np.radians(60.0))
        res = ev.evaluate(pos, np.ones(N), beta)
        assert res["peak_theta_deg"] == pytest.approx(60.0, abs=1.0)


class TestEvaluatorNonUniform:
    def test_nonuniform_runs_and_is_sane(self):
        ev = ArrayEvaluator(make_args())
        pos = np.array([0, 0.4, 0.9, 1.5, 2.2, 3.0, 3.9, 4.9])
        res = ev.evaluate(pos)
        assert res["N"] == 8
        assert 0 < res["directivity"] < 16
        assert res["peak_sidelobe_db"] < 0


class TestEvaluatorCSV:
    def test_read_geometry_csv_defaults(self, tmp_path):
        # amplitude/phase columns optional; should default to 1.0 / 0.0.
        p = tmp_path / "geom.csv"
        p.write_text("position_lambda\n0.0\n0.5\n1.0\n")
        pos, amps, ph = ArrayEvaluator.read_geometry_csv(str(p))
        assert np.allclose(pos, [0.0, 0.5, 1.0])
        assert np.allclose(amps, 1.0)
        assert np.allclose(ph, 0.0)

    def test_read_geometry_csv_full(self, tmp_path):
        p = tmp_path / "g.csv"
        p.write_text("position_lambda,amplitude,phase_deg\n0.0,1.0,0\n0.5,0.5,90\n")
        pos, amps, ph = ArrayEvaluator.read_geometry_csv(str(p))
        assert np.allclose(amps, [1.0, 0.5])
        assert np.allclose(ph, [0.0, np.pi / 2])

    def test_cli_evaluate(self, tmp_path):
        p = tmp_path / "g.csv"
        p.write_text("position_lambda,amplitude,phase_deg\n"
                     + "\n".join("%g,1.0,0" % (i * 0.5) for i in range(8)) + "\n")
        shell = ArrayCalculator(["evaluate", "-g", str(p), "--variable_return"])
        shell.main(shell.getArgs())
        res = shell.getCalcedParams()
        assert res["N"] == 8
        assert res["peak_theta_deg"] == pytest.approx(90.0, abs=0.5)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
