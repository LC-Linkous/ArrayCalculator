#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  triangular_array.py
#
#   Triangular linear array: amplitudes ramp linearly from the edges up
#   to the center. The distribution is the convolution of two uniform
#   arrays, which puts its first sidelobe at about -26 dB -- exactly the
#   square (in the pattern domain) of the uniform array's. A simple,
#   robust low-sidelobe taper. Reference: Balanis, Antenna Theory.
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees
import numpy as np

from array_common import ArrayCommon


class TriangularArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        # Linear taper, peak at center, tapering toward the edges. Using
        # (m+1) in the denominator keeps the edge elements nonzero (a true
        # zero-amplitude edge element is just an (N-2)-element array).
        m = (N - 1) / 2.0
        n = np.arange(N) - m
        return 1.0 - np.abs(n) / (m + 1.0)

    def hpbw(self, N):
        # Computed directly from the array factor (-3 dB crossings), the
        # same approach used by the uniform and Taylor arrays for accuracy
        # across N.
        amps = self.amplitudes(N)
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
        # Aperture (taper) efficiency times the uniform directivity:
        # D0 = (sum a)^2 / (N sum a^2) * 2 N d/lambda.
        amps = np.asarray(amps, dtype=float)
        eff = (amps.sum() ** 2) / (len(amps) * np.sum(amps ** 2))
        return eff * 2.0 * N * d_over_lambda

    # --- driver --------------------------------------------------------
    def triangular_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing

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
            self.plot_pattern(amps_norm, title="Triangular Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db