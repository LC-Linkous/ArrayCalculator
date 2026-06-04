#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  villeneuve_array.py
#
#   Villeneuve n-bar linear array synthesis. The discrete-array analog of
#   the Taylor n-bar line source: it holds the first nbar sidelobes near a
#   design level and lets the rest decay, but derives the pattern zeros
#   for a finite set of N elements directly (a Dolph-Tschebyscheff-style
#   discrete pattern), rather than sampling a continuous line-source
#   distribution as Taylor does. For small/moderate N this places the
#   near-in sidelobes more accurately than sampled Taylor; as N grows the
#   two converge. Reference: Villeneuve (1984); Balanis, Antenna Theory.
#
#   Construction: the discrete (N-element) array has N-1 pattern zeros.
#   Villeneuve keeps the inner zeros at the Dolph-Tschebyscheff (equal-
#   ripple) locations, stretched by a factor sigma so that the nbar-th
#   zero lands on the uniform-array integer location; the outer zeros are
#   left at the uniform-array integers n = nbar .. (so the far sidelobes
#   decay like a uniform array). The excitations are the inverse DFT of
#   the resulting zero pattern.
##--------------------------------------------------------------------\

from math import log10, sqrt, pi, radians, degrees
import numpy as np

from array_common import ArrayCommon


class VilleneuveArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def voltage_ratio(self, sll_db):
        return 10.0 ** (sll_db / 20.0)

    def amplitudes(self, N, sll_db, nbar):
        # An N-element array factor is a degree-(N-1) polynomial in
        # z = exp(j psi); it has exactly N-1 zeros on the unit circle,
        # symmetric about psi = pi. We place them and recover the element
        # excitations as the polynomial coefficients.
        R = self.voltage_ratio(sll_db)
        P = N - 1                      # number of pattern zeros

        # Uniform-array zeros are evenly spaced: psi_n = 2 pi n / N,
        # n = 1 .. N-1  (i.e. all the grating-free nulls of a uniform array).
        n = np.arange(1, N)            # 1 .. N-1
        psi_uniform = 2.0 * np.pi * n / N

        # Dolph-Tschebyscheff (equal-ripple) zeros for the same N: the nulls
        # of T_P(x0 cos(psi/2)). x0 cos(psi_m/2) = cos((2m-1) pi / (2P)).
        x0 = np.cosh(np.arccosh(R) / P)
        m = np.arange(1, P + 1)
        psi_dolph = 2.0 * np.arccos(
            np.clip(np.cos((2 * m - 1) * np.pi / (2 * P)) / x0, -1.0, 1.0))
        psi_dolph = np.sort(psi_dolph)     # ascending in (0, 2pi)

        # Villeneuve stretch sigma: scale the inner Dolph zeros so the nbar-th
        # one lands on the uniform integer location; keep uniform zeros beyond
        # nbar so the far sidelobes decay like a uniform array.
        sigma = psi_uniform[nbar - 1] / psi_dolph[nbar - 1]

        psi_zeros = psi_uniform.copy()
        # inner zeros 1..nbar-1 (and their mirror images about pi) take the
        # stretched-Dolph locations; the rest stay uniform.
        for k in range(nbar - 1):
            psi_zeros[k] = sigma * psi_dolph[k]
            psi_zeros[N - 2 - k] = 2.0 * np.pi - sigma * psi_dolph[k]

        zeros = np.exp(1j * psi_zeros)
        coeffs = np.poly(zeros)            # length N
        amps = np.abs(np.real(coeffs))
        amps = 0.5 * (amps + amps[::-1])   # enforce symmetry numerically
        return amps, R, sigma

    def directivity(self, N, d_over_lambda, amps):
        amps = np.asarray(amps, dtype=float)
        eff = (amps.sum() ** 2) / (len(amps) * np.sum(amps ** 2))
        return eff * 2.0 * N * d_over_lambda

    def hpbw(self, N, amps):
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
    def villeneuve_array_calculator(self):
        N = self.args.elements
        sll_db = self.args.sidelobe_level
        nbar = self.args.nbar
        d_over_lambda = self.args.spacing

        amps, R, sigma = self.amplitudes(N, sll_db, nbar)
        amps_norm = self.normalize(amps)

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
                              title="Villeneuve Array (N=%d, %g dB, nbar=%d)" % (N, sll_db, nbar),
                              sll_db=sll_db)

        if self.args.variable_return:
            return amps_norm, degrees(hpbw_rad), D0_db
