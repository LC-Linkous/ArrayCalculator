#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_normalization.py
#
#   Tests for the shared ArrayCommon.normalize() helper and the universal
#   --norm flag. The key invariants:
#     - normalization changes only the overall scale, never the ordering
#       or shape of the amplitude vector;
#     - 'edge' and 'center' results are scalar multiples of each other;
#     - the radiation pattern (and therefore the sidelobe level) is
#       identical regardless of --norm -- it is a cosmetic choice;
#     - the edge mode is guarded against near-zero edge samples (Blackman).
##--------------------------------------------------------------------\

import numpy as np
import pytest

from uniform_array import UniformArray
from binomial_array import BinomialArray
from triangular_array import TriangularArray
from blackman_array import BlackmanArray
from dolph_tschebyscheff import DolphTschebyscheff
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db


# A center-peaked vector to exercise the helper directly.
_PEAKED = np.array([1.042, 1.413, 2.051, 2.584, 2.887,
                    2.887, 2.584, 2.051, 1.413, 1.042])


class TestNormalizeHelper:
    def test_center_peak_is_one(self):
        u = UniformArray(make_args(elements=10))
        out = u.normalize(_PEAKED, mode="center")
        assert out.max() == pytest.approx(1.0)
        assert np.argmax(out) == np.argmax(_PEAKED)

    def test_edge_min_is_one(self):
        u = UniformArray(make_args(elements=10))
        out = u.normalize(_PEAKED, mode="edge")
        assert out.min() == pytest.approx(1.0)
        # edge-normalized Dolph-like array rises toward the center
        assert out[len(out) // 2] > out[0]

    def test_ordering_preserved(self):
        u = UniformArray(make_args(elements=10))
        e = u.normalize(_PEAKED, mode="edge")
        c = u.normalize(_PEAKED, mode="center")
        assert np.argsort(e).tolist() == np.argsort(_PEAKED).tolist()
        assert np.argsort(c).tolist() == np.argsort(_PEAKED).tolist()

    def test_edge_and_center_are_scalar_multiples(self):
        u = UniformArray(make_args(elements=10))
        e = u.normalize(_PEAKED, mode="edge")
        c = u.normalize(_PEAKED, mode="center")
        ratio = e / c
        assert np.allclose(ratio, ratio[0])

    def test_edge_guard_against_near_zero(self):
        # Blackman with the /N convention has a tiny (but nonzero) edge
        # sample; the guard must keep edge-normalization finite and
        # well-behaved rather than dividing by ~0.
        obj = BlackmanArray(make_args(elements=16))
        amps = obj.amplitudes(16)
        out = obj.normalize(amps, mode="edge")
        assert np.all(np.isfinite(out))
        assert np.argmax(out) == np.argmax(amps)

    def test_mode_defaults_to_args(self):
        # When mode is omitted, the helper reads args.norm.
        u_edge = UniformArray(make_args(elements=10, norm="edge"))
        u_ctr = UniformArray(make_args(elements=10, norm="center"))
        assert u_edge.normalize(_PEAKED).min() == pytest.approx(1.0)
        assert u_ctr.normalize(_PEAKED).max() == pytest.approx(1.0)


class TestNormalizationPatternInvariance:
    # The pattern must not depend on the normalization convention. If this
    # ever fails, normalization has stopped being a pure scalar scaling.
    @pytest.mark.parametrize("cls,kw", [
        (TriangularArray, {}),
        (BlackmanArray, {}),
        (BinomialArray, {}),
    ])
    def test_sidelobe_level_invariant(self, cls, kw):
        obj = cls(make_args(elements=20, **kw))
        amps = obj.amplitudes(20)
        sll_edge = peak_sidelobe_db(obj, obj.normalize(amps, mode="edge"))
        sll_ctr = peak_sidelobe_db(obj, obj.normalize(amps, mode="center"))
        assert sll_edge == pytest.approx(sll_ctr, abs=1e-6)

    def test_dolph_sidelobe_invariant(self):
        dt = DolphTschebyscheff(make_args(elements=12, sidelobe_level=25))
        amps, _, _ = dt.amplitudes(12, 25)
        e = peak_sidelobe_db(dt, dt.normalize(amps, mode="edge"))
        c = peak_sidelobe_db(dt, dt.normalize(amps, mode="center"))
        assert e == pytest.approx(c, abs=1e-6)
        assert e == pytest.approx(-25, abs=0.1)


class TestNormCLI:
    def test_norm_is_universal_flag(self):
        # --norm is accepted by a closed-form taper, not just Dolph.
        shell = ArrayCalculator(["triangular_array", "-N", "10",
                                 "--norm", "edge", "--variable_return"])
        shell.main(shell.getArgs())
        amps, _, _ = shell.getCalcedParams()
        assert amps.min() == pytest.approx(1.0)

    def test_taper_default_is_center(self):
        shell = ArrayCalculator(["triangular_array", "-N", "10",
                                 "--variable_return"])
        assert shell.getArgs().norm == "center"
        shell.main(shell.getArgs())
        amps, _, _ = shell.getCalcedParams()
        assert amps.max() == pytest.approx(1.0)

    def test_dolph_default_is_edge(self):
        shell = ArrayCalculator(["dolph_tschebyscheff", "-N", "10",
                                 "-sll", "26", "--variable_return"])
        assert shell.getArgs().norm == "edge"
        shell.main(shell.getArgs())
        amps, _, _, _ = shell.getCalcedParams()
        assert amps.min() == pytest.approx(1.0)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))