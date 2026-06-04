#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  tests/test_helpers.py
#
#   Shared test helpers ONLY. This module is imported by the per-array
#   test files (test_uniform_array.py, test_binomial_array.py, ...) and
#   by test_array_calculator.py for the cross-cutting suites. It defines
#   no test classes itself, so pytest collects nothing from it -- despite
#   the test_ prefix, which is kept so the import path stays simple.
#
#   If you add a helper used by more than one test file, it belongs here.
##--------------------------------------------------------------------\

from types import SimpleNamespace

import numpy as np


def make_args(**kwargs):
    """Build an args namespace with the defaults the array classes expect.

    Mirrors the production CLI defaults from array_calculator._add_common_args
    (notably norm='center'); pass norm='edge' explicitly for Dolph-style
    edge-normalized checks. Override any field via keyword.
    """
    defaults = dict(
        elements=None, frequency=None, spacing=0.5, scan=90.0,
        unit="centimeter", csv=None, plot=False, verbose=False,
        variable_return=False, sidelobe_level=None, norm="center",
        plot_style="both", save=None,
    )
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


def peak_sidelobe_db(obj, amps):
    """Highest sidelobe level (dB) of the array's actual radiation pattern.

    Finds every interior local maximum of the normalized pattern below the
    main beam (< -0.5 dB) and returns the largest. This is the model-
    independent ground truth used throughout the suite to check that a
    synthesis method's sidelobes land where the theory says they should.
    """
    theta, af_norm, af_db = obj.pattern_sweep(amps, n_points=20001)
    peaks = [af_db[i] for i in range(1, len(af_db) - 1)
             if af_db[i] >= af_db[i - 1] and af_db[i] >= af_db[i + 1]
             and af_db[i] < -0.5]
    return max(peaks) if peaks else -np.inf