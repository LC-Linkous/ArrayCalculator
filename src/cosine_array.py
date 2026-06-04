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
            self.value_print("HPBW", degrees(hpbw_rad), "deg")
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm, title="Cosine Array (N=%d)" % N,
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db