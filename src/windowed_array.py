#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  windowed_array.py
#
#   STAGED / NOT YET WIRED IN.
#
#   This module is the future consolidation of the closed-form amplitude
#   tapers that currently live in their own one-file-per-array modules
#   (triangular_array.py, cosine_array.py, cosine_squared_array.py,
#   hann_array.py, hamming_array.py, blackman_array.py). Every one of
#   those tapers differs ONLY in its amplitude formula; the HPBW (from
#   the array factor) and directivity (aperture efficiency) machinery is
#   identical. This class collapses all six into a single taper table so
#   the shared logic lives in one place.
#
#   It is intentionally left UNREGISTERED in array_calculator.py for now.
#   The standalone classes remain the source of truth and the ones the
#   tests import. When we are ready to consolidate, point the CLI here
#   and retire the six modules. Nothing imports this yet, so it has no
#   effect on the current calculator.
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees, pi
import numpy as np

from array_common import ArrayCommon


# Taper formulas, keyed by name. Each maps N -> raw amplitude array.
# Kept here as the single definition of each closed-form distribution.
def _triangular(N):
    m = (N - 1) / 2.0
    n = np.arange(N) - m
    return 1.0 - np.abs(n) / (m + 1.0)


def _cosine(N):
    n = np.arange(N) - (N - 1) / 2.0
    return np.abs(np.cos(pi * n / N))


def _cosine_squared(N):
    n = np.arange(N) - (N - 1) / 2.0
    return np.cos(pi * n / N) ** 2


def _hann(N):
    x = (np.arange(N) - (N - 1) / 2.0) / N
    return 0.5 + 0.5 * np.cos(2.0 * pi * x)


def _hamming(N):
    x = (np.arange(N) - (N - 1) / 2.0) / N
    return 0.54 + 0.46 * np.cos(2.0 * pi * x)


def _blackman(N):
    x = (np.arange(N) - (N - 1) / 2.0) / N
    return 0.42 + 0.5 * np.cos(2.0 * pi * x) + 0.08 * np.cos(4.0 * pi * x)


TAPERS = {
    "triangular": ("Triangular", _triangular),
    "cosine": ("Cosine", _cosine),
    "cosine_squared": ("Cosine-Squared", _cosine_squared),
    "hann": ("Hann", _hann),
    "hamming": ("Hamming", _hamming),
    "blackman": ("Blackman", _blackman),
}


class WindowedArray(ArrayCommon):
    def __init__(self, args, taper=None):
        super().__init__(args)
        self.args = args
        # taper may be passed explicitly or read from args.taper later,
        # once this module is wired into the CLI.
        self.taper = taper or getattr(args, "taper", None)

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N, taper=None):
        key = taper or self.taper
        if key not in TAPERS:
            raise ValueError("unknown taper %r; choose from %s"
                             % (key, ", ".join(sorted(TAPERS))))
        return np.abs(TAPERS[key][1](N))

    def hpbw(self, N, taper=None):
        amps = self.amplitudes(N, taper)
        theta, af_norm, _ = self.pattern_sweep(amps, n_points=40001)
        half = 1.0 / sqrt(2.0)
        peak = np.argmax(af_norm)
        left = peak
        while left > 0 and af_norm[left] > half:
            left -= 1
        right = peak
        while right < len(af_norm) - 1 and af_norm[right] > half:
            right += 1
        return radians(theta[right] - theta[left])

    def directivity(self, N, d_over_lambda, amps):
        amps = np.asarray(amps, dtype=float)
        eff = (amps.sum() ** 2) / (len(amps) * np.sum(amps ** 2))
        return eff * 2.0 * N * d_over_lambda

    # --- driver --------------------------------------------------------
    def windowed_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing
        label = TAPERS[self.taper][0]

        amps = self.amplitudes(N)
        amps_norm = self.normalize(amps)
        hpbw_rad = self.hpbw(N)
        D0 = self.directivity(N, d_over_lambda, amps)
        D0_db = 10.0 * log10(D0)

        if not self.args.variable_return:
            self.value_print("N", N)
            self.list_print("Amplitudes (normalized)", amps_norm)
            self.report_spacing(d_over_lambda)
            self.value_print("HPBW", degrees(hpbw_rad), "deg")
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm,
                              title="%s Array (N=%d)" % (label, N),
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db