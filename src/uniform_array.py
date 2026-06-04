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
        # Broadside HPBW for a uniform array (Balanis). For a large array
        # this is well approximated from the first-null/half-power relation;
        # we compute it directly from the array factor for accuracy across N.
        amps = self.amplitudes(N)
        theta, af_norm, _ = self.pattern_sweep(amps, n_points=40001)
        # find the two half-power (-3 dB, i.e. 1/sqrt(2)) crossings around the
        # broadside peak at 90 deg
        half = 1.0 / sqrt(2.0)
        peak_idx = np.argmax(af_norm)
        # walk left and right to the half-power points
        left = peak_idx
        while left > 0 and af_norm[left] > half:
            left -= 1
        right = peak_idx
        while right < len(af_norm) - 1 and af_norm[right] > half:
            right += 1
        return radians(theta[right] - theta[left])

    def directivity(self, N, d_over_lambda):
        # D0 = 2 N d / lambda for a broadside uniform array (large-N form);
        # for d = lambda/2 this is simply N.
        return 2.0 * N * d_over_lambda

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
            self.value_print("HPBW", degrees(hpbw_rad), "deg")
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Uniform Array (N=%d)" % N)

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db