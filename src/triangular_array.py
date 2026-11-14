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
        return self.hpbw_rad(amps)

    def directivity(self, N, d_over_lambda, amps):
        # Exact pattern-integral directivity: correct for any spacing
        # (including grating lobes) and any --scan angle.
        return self.directivity_integral(amps, d_over_lambda=d_over_lambda)

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
            self.print_hpbw(hpbw_rad)
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")
            self.warn_grating_lobes(amps)

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Triangular Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, self.hpbw_out(hpbw_rad), D0_db