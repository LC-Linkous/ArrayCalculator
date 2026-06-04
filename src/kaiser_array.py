#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  kaiser_array.py
#
#   Kaiser (Kaiser-Bessel) linear array: a one-parameter taper whose
#   shape parameter beta trades main-beam width against sidelobe level
#   continuously. It is the closed-form taper that comes closest to the
#   prolate-spheroidal optimum, and the natural bridge between the fixed
#   window tapers (triangular, Hann, ...) and the parametric synthesis
#   methods (Dolph, Taylor): small beta -> near-uniform (high sidelobes,
#   narrow beam); large beta -> heavy taper (low sidelobes, wide beam).
#
#       a(n) = I0( beta * sqrt(1 - (2 n / (N-1))^2) ) / I0(beta)
#
#   with n = -(N-1)/2 .. (N-1)/2 and I0 the modified Bessel function of
#   the first kind, order 0. Reference: Balanis, Antenna Theory; Kaiser
#   (1974); Harris (1978).
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees
import numpy as np

from array_common import ArrayCommon


def _i0(x):
    # Modified Bessel I0 via the standard series; vectorized over x. Converges
    # quickly for the beta range of interest (0..~12). Avoids a SciPy dep.
    x = np.asarray(x, dtype=float)
    total = np.ones_like(x)
    term = np.ones_like(x)
    k = 1
    while True:
        term = term * (x / (2.0 * k)) ** 2
        total = total + term
        if np.all(term < 1e-12 * np.maximum(total, 1e-12)):
            break
        k += 1
        if k > 200:
            break
    return total


class KaiserArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N, beta=None):
        if beta is None:
            beta = getattr(self.args, "beta", 6.0)
        n = np.arange(N) - (N - 1) / 2.0
        ratio = 2.0 * n / (N - 1) if N > 1 else np.zeros(N)
        return _i0(beta * np.sqrt(np.clip(1.0 - ratio ** 2, 0.0, None))) / _i0(beta)

    def hpbw(self, N, beta=None):
        amps = self.amplitudes(N, beta)
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
    def kaiser_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing
        beta = getattr(self.args, "beta", 6.0)

        amps = self.amplitudes(N, beta)
        amps_norm = self.normalize(amps)
        hpbw_rad = self.hpbw(N, beta)
        D0 = self.directivity(N, d_over_lambda, amps)
        D0_db = 10.0 * log10(D0)

        if not self.args.variable_return:
            self.value_print("N", N)
            self.value_print("beta", beta)
            self.list_print("Amplitudes (normalized)", amps_norm)
            self.report_spacing(d_over_lambda)
            self.value_print("HPBW", degrees(hpbw_rad), "deg")
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Kaiser Array (N=%d, beta=%g)" % (N, beta),
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db
