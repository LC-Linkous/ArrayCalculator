#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_window_def.py
#
#   --window-def selects the sampling of the Hann / Hamming / Blackman
#   tapers: 'aperture' x = n/N (default, unchanged behavior) or
#   'symmetric' x = n/(N-1) (DSP-style window, edges on the endpoints).
##--------------------------------------------------------------------\

import numpy as np
import pytest

from hann_array import HannArray
from hamming_array import HammingArray
from blackman_array import BlackmanArray
from cosine_squared_array import CosineSquaredArray
from array_calculator import ArrayCalculator
from test_helpers import make_args, peak_sidelobe_db

N = 16


@pytest.mark.parametrize("cls", [HannArray, HammingArray, BlackmanArray])
def test_default_is_aperture_and_unchanged(cls):
    a = cls(make_args(elements=N)).amplitudes(N)
    b = cls(make_args(elements=N, window_def="aperture")).amplitudes(N)
    np.testing.assert_allclose(a, b)
    x = (np.arange(N) - (N - 1) / 2) / N
    if cls is HannArray:
        np.testing.assert_allclose(a, 0.5 + 0.5 * np.cos(2 * np.pi * x))


@pytest.mark.parametrize("cls,edge", [(HannArray, 0.0), (HammingArray, 0.08),
                                      (BlackmanArray, 0.0)])
def test_symmetric_hits_window_endpoints(cls, edge):
    a = cls(make_args(elements=N, window_def="symmetric")).amplitudes(N)
    assert a[0] == pytest.approx(edge, abs=1e-12)
    assert a[-1] == pytest.approx(edge, abs=1e-12)
    np.testing.assert_allclose(a, a[::-1])


def test_hann_equals_cos2_only_in_aperture_mode():
    c2 = CosineSquaredArray(make_args(elements=N)).amplitudes(N)
    ap = HannArray(make_args(elements=N)).amplitudes(N)
    sy = HannArray(make_args(elements=N, window_def="symmetric")).amplitudes(N)
    np.testing.assert_allclose(ap, c2, atol=1e-12)
    assert not np.allclose(sy, c2)


def test_symmetric_hann_lower_sidelobes_wider_beam():
    ap = HannArray(make_args(elements=N))
    sy = HannArray(make_args(elements=N, window_def="symmetric"))
    a_ap, a_sy = ap.amplitudes(N), sy.amplitudes(N)
    assert peak_sidelobe_db(sy, a_sy) < peak_sidelobe_db(ap, a_ap)
    assert sy.hpbw_deg(a_sy) > ap.hpbw_deg(a_ap)


def test_cli_flag():
    shell = ArrayCalculator(["hann_array", "-N", "8", "--window-def", "symmetric",
                             "--norm", "edge", "--variable_return"])
    shell.main(shell.getArgs())
    amps, _, _ = shell.getCalcedParams()
    assert np.all(np.isfinite(amps))          # edge-norm skips the zero edges
    with pytest.raises(SystemExit):
        ArrayCalculator(["cosine_array", "-N", "8", "--window-def", "symmetric"])
