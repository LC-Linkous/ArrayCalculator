#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  taylor_array.py
#
#   Taylor (n-bar) linear array synthesis. Holds the near-in sidelobes
#   close to a specified level while letting the far-out sidelobes decay
#   (unlike Dolph-Tschebyscheff, whose sidelobes are all equal). This is
#   the distribution most often used in practice for large arrays because
#   it trades a slightly wider beam for better aperture efficiency and
#   more realistic, decaying sidelobes. Reference: Balanis, Antenna Theory.
#
#   The parameter nbar sets how many near-in sidelobes are held near the
#   design level before the natural (uniform-like) decay takes over.
##--------------------------------------------------------------------\

from math import log10, sqrt, pi, radians, degrees
import numpy as np

from array_common import ArrayCommon


class TaylorArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def voltage_ratio(self, sll_db):
        return 10.0 ** (sll_db / 20.0)

    def _zero(self, n, nbar, A, sigma):
        # nth zero of the Taylor pattern: scaled for n < nbar, integer beyond.
        if n < nbar:
            return sigma * sqrt(A ** 2 + (n - 0.5) ** 2)
        return float(n)

    def amplitudes(self, N, sll_db, nbar):
        # Taylor line-source aperture distribution sampled at the N element
        # positions. The distribution is
        #   I(p) = 1 + 2 * sum_{m=1}^{nbar-1} F(m) cos(2 pi m p)
        # where F(m) are the Taylor coefficients built from the pattern zeros
        # and p is the normalized element position. Verified numerically to
        # place the near-in sidelobes at the design level with decaying tails.
        R = self.voltage_ratio(sll_db)
        A = np.arccosh(R) / pi
        sigma = nbar / sqrt(A ** 2 + (nbar - 0.5) ** 2)

        Fm = []
        for m in range(1, nbar):
            num = 1.0
            for n in range(1, nbar):
                num *= (1.0 - (m ** 2) / (self._zero(n, nbar, A, sigma) ** 2))
            den = 1.0
            for n in range(1, nbar):
                if n != m:
                    den *= (1.0 - (m ** 2) / (n ** 2))
            Fm.append(((-1) ** (m + 1)) * 0.5 * num / den)
        Fm = np.array(Fm)

        n_idx = np.arange(N) - (N - 1) / 2.0
        p = n_idx / N
        amps = np.ones(N)
        for m in range(1, nbar):
            amps += 2.0 * Fm[m - 1] * np.cos(2.0 * pi * m * p)
        return np.abs(amps), R, sigma

    def directivity(self, N, d_over_lambda, amps):
        # Aperture (taper) efficiency from the excitation, times the uniform
        # directivity 2 N d/lambda. D0 = (sum a)^2 / (N sum a^2) * 2 N d/lambda.
        amps = np.asarray(amps, dtype=float)
        eff = (amps.sum() ** 2) / (len(amps) * np.sum(amps ** 2))
        return eff * 2.0 * N * d_over_lambda

    def hpbw(self, N, amps):
        # Computed directly from the array factor (-3 dB crossings).
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

    # --- driver --------------------------------------------------------
    def taylor_array_calculator(self):
        N = self.args.elements
        sll_db = self.args.sidelobe_level
        nbar = self.args.nbar
        d_over_lambda = self.args.spacing

        amps, R, sigma = self.amplitudes(N, sll_db, nbar)
        amps_norm = self.normalize(amps)  # Taylor conventionally center-normalized (the default)

        hpbw_rad = self.hpbw(N, amps)
        D0 = self.directivity(N, d_over_lambda, amps)
        D0_db = 10.0 * log10(D0)

        if not self.args.variable_return:
            self.value_print("N", N)
            self.value_print("Sidelobe level", sll_db, "dB")
            self.value_print("nbar", nbar)
            if self.args.verbose:
                self.value_print("R (voltage ratio)", R)
                self.value_print("sigma", sigma)
                self.list_print("Amplitudes (raw)", amps)
            self.list_print("Amplitudes (normalized)", amps_norm)
            self.report_spacing(d_over_lambda)
            self.value_print("HPBW", degrees(hpbw_rad), "deg")
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm,
                              title="Taylor Array (N=%d, %g dB, nbar=%d)" % (N, sll_db, nbar),
                              sll_db=sll_db)

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db