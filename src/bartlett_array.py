#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  bartlett_array.py
#
#   Bartlett linear array: the triangular (linear) taper in its
#   zero-endpoint form,
#
#       a(n) = 1 - |2 n / (N-1)|,   n = -(N-1)/2 .. (N-1)/2
#
#   so the two edge elements go to exactly zero. This is the textbook
#   "Bartlett window" and differs from triangular_array.py, which uses an
#   (m+1) denominator to keep the edge elements nonzero (a true zero-edge
#   element is electrically an (N-2)-element array). The two are the same
#   distribution apart from that endpoint convention; both place the first
#   sidelobe near -26 dB. Provided separately so the endpoint distinction
#   is explicit. Reference: Balanis, Antenna Theory; Harris (1978).
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees
import numpy as np

from array_common import ArrayCommon


class BartlettArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        n = np.arange(N) - (N - 1) / 2.0
        ratio = 2.0 * n / (N - 1) if N > 1 else np.zeros(N)
        return np.clip(1.0 - np.abs(ratio), 0.0, None)

    def hpbw(self, N):
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
        amps = np.asarray(amps, dtype=float)
        eff = (amps.sum() ** 2) / (len(amps) * np.sum(amps ** 2))
        return eff * 2.0 * N * d_over_lambda

    # --- driver --------------------------------------------------------
    def bartlett_array_calculator(self):
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
            self.plot_pattern(amps_norm, title="Bartlett Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db
