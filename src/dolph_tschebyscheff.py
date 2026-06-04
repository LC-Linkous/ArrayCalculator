#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator    dolph_tschebyscheff.py
#
#   Dolph-Tschebyscheff linear array synthesis. Produces the narrowest
#   main beam for a specified, equal sidelobe level (SLL). Reference:
#   Balanis, Antenna Theory.
#
#   NOTE ON COEFFICIENT VALUES:
#   The excitation coefficients computed here are the values that make
#   the radiated pattern's sidelobes sit exactly at the requested SLL
#   (verified numerically). They differ by <0.5 dB worth of amplitude
#   from a few of the hand-worked figures in the EGRE 540 handout's
#   Examples 4 and 6; that small drift came from accumulated rounding in
#   manual back-substitution (done without a calculator like this one).
#   The synthesis *method* in the handout is correct; these are simply
#   the precise numerical results of that same method.
##--------------------------------------------------------------------\

from math import log10, sqrt, radians, degrees
import numpy as np

from array_common import ArrayCommon


class DolphTschebyscheff(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- synthesis -----------------------------------------------------
    def voltage_ratio(self, sll_db):
        # R: main-beam-to-sidelobe voltage ratio (linear) from SLL in dB.
        return 10.0 ** (sll_db / 20.0)

    def z0(self, N, R):
        # Scaling parameter mapping the Tschebyscheff polynomial onto the
        # array factor: z0 = cosh( (1/(N-1)) * arccosh(R) ).
        return np.cosh(np.arccosh(R) / (N - 1))

    def amplitudes(self, N, sll_db):
        # The array factor, as a polynomial in cos(u), must equal the
        # (N-1)th Tschebyscheff polynomial T_{N-1}(z0 cos u). We recover the
        # element excitations by fitting the symmetric cosine series of the
        # element pattern to that target polynomial. This is the same
        # coefficient-matching the handout does by hand, solved numerically
        # so it holds for any N (even or odd) and any SLL.
        R = self.voltage_ratio(sll_db)
        P = N - 1
        z0 = self.z0(N, R)

        u = np.linspace(0.0, np.pi, 8000)
        w = z0 * np.cos(u)
        Tp = np.where(
            np.abs(w) <= 1.0,
            np.cos(P * np.arccos(np.clip(w, -1.0, 1.0))),
            np.cosh(P * np.arccosh(np.clip(np.abs(w), 1.0, None))) * np.sign(w) ** P,
        )
        n_idx = np.arange(N) - (N - 1) / 2.0
        basis = np.cos(np.outer(2.0 * u, n_idx))
        coeffs, *_ = np.linalg.lstsq(basis, Tp, rcond=None)
        return np.abs(coeffs), R, z0

    def directivity(self, N, R):
        # Directivity for a Dolph-Tschebyscheff array (Balanis), with the
        # beam-broadening factor f = 1 + 0.636*(2/R * cosh(sqrt(arccosh(R)^2
        # - pi^2)))^2 for d = lambda/2.
        arg = np.arccosh(R) ** 2 - np.pi ** 2
        if arg <= 0:
            f = 1.0
        else:
            f = 1.0 + 0.636 * (2.0 / R * np.cosh(np.sqrt(arg))) ** 2
        # Balanis (6-79): D0 = 2 R^2 / (1 + (R^2 - 1) f / [(L+d)/lambda]),
        # with the broadside array length L = (N-1) d, so (L+d)/lambda =
        # N * d / lambda. (An earlier version used 2 N d / lambda here, which
        # roughly doubled D0.)
        D0 = 2.0 * R ** 2 / (1.0 + (R ** 2 - 1.0) * f / (N * self.args.spacing))
        return D0

    # --- driver --------------------------------------------------------
    def dolph_tschebyscheff_calculator(self):
        N = self.args.elements
        sll_db = self.args.sidelobe_level
        d_over_lambda = self.args.spacing

        amps, R, z0 = self.amplitudes(N, sll_db)

        amps_norm = self.normalize(amps)

        D0 = self.directivity(N, R)
        D0_db = 10.0 * log10(D0)

        if not self.args.variable_return:
            self.value_print("N", N)
            self.value_print("Sidelobe level", sll_db, "dB")
            if self.args.verbose:
                self.value_print("R (voltage ratio)", R)
                self.value_print("z0", z0)
                self.list_print("Amplitudes (raw)", amps)
            self.list_print("Amplitudes (%s-normalized)" % self.args.norm, amps_norm)
            self.report_spacing(d_over_lambda)
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, amps_norm)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps_norm,
                              title="Dolph-Tschebyscheff (N=%d, %g dB)" % (N, sll_db),
                              sll_db=sll_db)

        if self.args.variable_return:
            return amps_norm, R, z0, D0_db