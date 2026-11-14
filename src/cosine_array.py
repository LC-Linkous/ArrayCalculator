#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  cosine_array.py
#
#   Cosine (half-cosine) linear array: amplitudes follow a single cosine
#   lobe across the aperture, zero-ish at the edges and peak at center.
#   First sidelobe lands at about -23 dB -- a modest, smooth taper that
#   trades a little beamwidth for a ~10 dB sidelobe improvement over a
#   uniform array. Reference: Balanis, Antenna Theory.
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees, pi
import numpy as np

from array_common import ArrayCommon


class CosineArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N):
        # Single cosine lobe sampled at the element positions: a(n) =
        # cos(pi * p), p in (-1/2, 1/2) across the aperture.
        n = np.arange(N) - (N - 1) / 2.0
        return np.abs(np.cos(pi * n / N))

    def hpbw(self, N):
        amps = self.amplitudes(N)
        return self.hpbw_rad(amps)

    def directivity(self, N, d_over_lambda, amps):
        # Exact pattern-integral directivity: correct for any spacing
        # (including grating lobes) and any --scan angle.
        return self.directivity_integral(amps, d_over_lambda=d_over_lambda)

    # --- driver --------------------------------------------------------
    def cosine_array_calculator(self):
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
            self.plot_pattern(amps_norm, title="Cosine Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, self.hpbw_out(hpbw_rad), D0_db