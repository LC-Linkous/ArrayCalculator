#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  hamming_array.py
#
#   Hamming linear array: a(n) = 0.54 + 0.46 cos(2 pi p). The pedestal is
#   tuned to cancel the nearest sidelobe, pushing the peak sidelobe down
#   to roughly -40 to -42 dB (it improves with N). Unlike Hann, the far
#   sidelobes do not decay as fast because the edge amplitudes are not
#   driven to zero. Reference: Balanis, Antenna Theory; Harris (1978).
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees, pi
import numpy as np

from array_common import ArrayCommon


class HammingArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        # Raised cosine with the classic 0.54/0.46 Hamming pedestal.
        x = self.window_x(N)   # n/N (default) or n/(N-1); see --window-def
        return 0.54 + 0.46 * np.cos(2.0 * pi * x)

    def hpbw(self, N):
        amps = self.amplitudes(N)
        return self.hpbw_rad(amps)

    def directivity(self, N, d_over_lambda, amps):
        # Exact pattern-integral directivity: correct for any spacing
        # (including grating lobes) and any --scan angle.
        return self.directivity_integral(amps, d_over_lambda=d_over_lambda)

    # --- driver --------------------------------------------------------
    def hamming_array_calculator(self):
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
            self.plot_pattern(amps_norm, title="Hamming Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, self.hpbw_out(hpbw_rad), D0_db