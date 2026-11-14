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
        return self.hpbw_rad(amps)

    def directivity(self, N, d_over_lambda, amps):
        # Exact pattern-integral directivity: correct for any spacing
        # (including grating lobes) and any --scan angle.
        return self.directivity_integral(amps, d_over_lambda=d_over_lambda)

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
            self.print_hpbw(hpbw_rad)
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")
            self.warn_grating_lobes(amps)

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Bartlett Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, self.hpbw_out(hpbw_rad), D0_db
