#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator    binomial_array.py
#
#   Binomial linear array synthesis. Amplitudes follow the binomial
#   coefficients (Pascal's triangle), giving a pattern with no sidelobes
#   at the cost of a wider main beam. Reference: Balanis, Antenna Theory.
##--------------------------------------------------------------------\

from math import comb, sqrt, radians, degrees, log10, cos, pi
import numpy as np
from pint import UnitRegistry  # pip install pint
ureg = UnitRegistry()

from array_common import ArrayCommon


class BinomialArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        # Binomial coefficients C(N-1, n), i.e. the (N-1)th row of Pascal's
        # triangle. These are the element excitation amplitudes.
        return np.array([comb(N - 1, n) for n in range(N)], dtype=float)

    def hpbw(self, N):
        # Half-power beamwidth for a binomial array with d = lambda/2,
        # HPBW ~= 1.06 / sqrt(N-1) radians (Balanis).
        return 1.06 / sqrt(N - 1)

    def directivity(self, N, d_over_lambda):
        # Balanis' closed-form *approximation* for the binomial array:
        # D0 ~= 1.77 * sqrt(1 + 2 L / lambda), L = (N-1) d. For d = lambda/2
        # this reduces to 1.77 * sqrt(N). It is a large-array fit and runs a
        # few percent high versus a direct pattern integral (e.g. N=6 gives
        # 4.34 vs ~4.06 integrated); kept here to match the EGRE 540 notes.
        L_over_lambda = (N - 1) * d_over_lambda
        return 1.77 * sqrt(1.0 + 2.0 * L_over_lambda)

    # --- driver --------------------------------------------------------
    def binomial_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing  # fraction of a wavelength

        amps = self.amplitudes(N)
        amps_norm = self.normalize(amps)
        hpbw_rad = self.hpbw(N)
        D0 = self.directivity(N, d_over_lambda)
        D0_db = 10.0 * log10(D0)

        if not self.args.variable_return:
            self.value_print("N", N)
            if self.args.verbose:
                self.list_print("Amplitudes (raw)", amps)
            self.list_print("Amplitudes (normalized)", amps_norm)
            self.report_spacing(d_over_lambda)        # prints lambda + d if -f given
            self.value_print("HPBW", degrees(hpbw_rad), "deg")
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        # pattern emission (csv now; plot is a future hook)
        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Binomial Array (N=%d)" % N)

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db