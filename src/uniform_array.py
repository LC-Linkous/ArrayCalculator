#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  uniform_array.py
#
#   Uniform linear array: all elements excited equally. This is the
#   reference case in the EGRE 540 notes -- narrowest beamwidth and
#   highest directivity of the classic distributions, at the cost of the
#   highest sidelobes (~ -13.2 dB). Reference: Balanis, Antenna Theory.
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees, asin, pi
import numpy as np

from array_common import ArrayCommon


class UniformArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        # Every element fed equally.
        return np.ones(N, dtype=float)

    def hpbw(self, N, d_over_lambda):
        # Read directly from the array factor (-3 dB crossings), so it holds
        # for any N, spacing and --scan angle, including endfire.
        return self.hpbw_rad(self.amplitudes(N))

    def directivity(self, N, d_over_lambda):
        # Exact pattern-integral directivity. Equals N at d = lambda/2
        # broadside and ~2 N d / lambda below the grating-lobe onset, but
        # stays correct when grating lobes appear or the beam is steered
        # (endfire roughly doubles it).
        return self.directivity_integral(self.amplitudes(N),
                                         d_over_lambda=d_over_lambda)

    # --- driver --------------------------------------------------------
    def uniform_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing

        amps = self.amplitudes(N)
        amps_norm = self.normalize(amps)
        hpbw_rad = self.hpbw(N, d_over_lambda)
        D0 = self.directivity(N, d_over_lambda)
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
            self.plot_pattern(amps_norm, title="Uniform Array (N=%d)" % N)

        if self.args.variable_return:
            return amps_norm, self.hpbw_out(hpbw_rad), D0_db