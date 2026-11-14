#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  cosine_squared_array.py
#
#   Cosine-squared (raised-cosine) linear array: the square of the cosine
#   taper. Smoother edges than plain cosine push the first sidelobe down
#   to about -31 dB at the cost of a wider main beam. Mathematically this
#   is the same distribution as the Hann window (cos^2 == 0.5 + 0.5 cos),
#   provided here under its antenna-theory name. Reference: Balanis,
#   Antenna Theory.
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees, pi
import numpy as np

from array_common import ArrayCommon


class CosineSquaredArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        # Square of the single cosine lobe.
        n = np.arange(N) - (N - 1) / 2.0
        return np.cos(pi * n / N) ** 2

    def hpbw(self, N):
        amps = self.amplitudes(N)
        return self.hpbw_rad(amps)

    def directivity(self, N, d_over_lambda, amps):
        # Exact pattern-integral directivity: correct for any spacing
        # (including grating lobes) and any --scan angle.
        return self.directivity_integral(amps, d_over_lambda=d_over_lambda)

    # --- driver --------------------------------------------------------
    def cosine_squared_array_calculator(self):
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
            self.print_hpbw(hpbw_rad)
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")
            self.warn_grating_lobes(amps)

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Cosine-Squared Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, self.hpbw_out(hpbw_rad), D0_db