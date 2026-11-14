#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  woodward_lawson_array.py
#
#   Woodward-Lawson shaped-beam synthesis. Unlike every other method in
#   this calculator (which produce a pencil beam with controlled
#   sidelobes), Woodward-Lawson designs an array to APPROXIMATE an
#   arbitrary desired pattern shape. It samples the desired pattern at a
#   set of orthogonal beam positions and sums one "composing" (sinc-like)
#   beam per sample, each steered to its sample point and weighted by the
#   desired amplitude there. The element excitations are the superposition
#   of those composing beams.
#
#   Two standard target shapes are built in:
#     * flat-top  : uniform amplitude across scan +/- sector, floor outside
#     * cosecant^2: classic one-sided radar ground-coverage pattern. Power
#                   ~ csc^2(elevation), peak at --scan, filling --sector
#                   degrees on the side away from broadside (see _target)
#
#   The excitations are complex (amplitude AND phase); every output
#   (pattern CSV, plot, directivity) uses the full complex excitation.
#
#   Reference: Woodward & Lawson (1948); Balanis, Antenna Theory.
##--------------------------------------------------------------------\

from math import log10, cos, pi, radians
import numpy as np

from array_common import ArrayCommon


class WoodwardLawsonArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- desired pattern targets --------------------------------------
    def _target(self, theta_deg):
        # Desired |pattern| (linear amplitude) as a function of theta (deg),
        # theta measured from the array axis.
        shape = getattr(self.args, "shape", "flat_top")
        center = getattr(self.args, "scan", 90.0)
        sector = getattr(self.args, "sector", 30.0)
        floor = getattr(self.args, "floor", 0.0)        # linear floor outside
        th = np.asarray(theta_deg, dtype=float)
        out = np.full_like(th, floor)
        if shape == "flat_top":
            # sector is the half-width around the scan angle
            inside = (th >= center - sector) & (th <= center + sector)
            out[inside] = 1.0
        elif shape == "cosecant_squared":
            # Treat broadside (theta = 90) as the horizon, so the elevation
            # of direction theta is eps = 90 - theta and sin(eps) = cos(theta).
            # The beam peaks at theta = scan (eps0 = 90 - scan, which must be
            # nonzero) and the amplitude falls as
            #   csc(eps) / csc(eps0) = cos(scan) / cos(theta)
            # (power ~ csc^2) across a fill region `sector` degrees wide on
            # the side away from the horizon: theta in [scan - sector, scan]
            # for scan < 90 (looking up), or [scan, scan + sector] for
            # scan > 90 (looking down). The opposite side of broadside gets
            # only the floor, so the target is one-sided by construction.
            if np.isclose(center, 90.0):
                raise ValueError("cosecant_squared needs --scan away from 90 "
                                 "(the peak must sit off the horizon, e.g. --scan 80)")
            if center < 90.0:
                lo, hi = max(center - sector, 0.0), center
            else:
                lo, hi = center, min(center + sector, 180.0)
            inside = (th >= lo) & (th <= hi)
            val = cos(radians(center)) / np.cos(np.radians(th[inside]))
            out[inside] = np.clip(val, 0.0, 1.0)
        else:
            raise ValueError("unknown shape %r" % shape)
        return out

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N, d_over_lambda=None):
        # Composing-beam sample points: an N-element array has N orthogonal
        # beams, spaced evenly in cos(theta). We sample the desired pattern
        # at those angles and superpose one steered composing beam per sample.
        # Returns COMPLEX excitations.
        if d_over_lambda is None:
            d_over_lambda = self.args.spacing
        m = np.arange(N) - (N - 1) / 2.0
        # orthogonal sample directions in cos(theta) space:
        #   cos(theta_m) = m / (N * d/lambda)
        cos_tm = m / (N * d_over_lambda)
        valid = np.abs(cos_tm) <= 1.0                  # only real (visible) angles
        theta_m = np.degrees(np.arccos(np.clip(cos_tm, -1.0, 1.0)))
        b_m = np.where(valid, self._target(theta_m), 0.0)   # desired samples
        if not np.any(b_m):
            raise ValueError(
                "no composing beam falls inside the shaped sector: the N=%d "
                "sample directions are %.1f deg apart near broadside; widen "
                "--sector, raise -N, or set a nonzero --floor" %
                (N, np.degrees(1.0 / (N * d_over_lambda))))

        # The array factor is AF = sum_n a_n exp(+j 2 pi x_n cos(theta)) with
        # x_n = n d/lambda, so a composing beam that peaks at cos(theta_m)
        # needs the conjugate phase:
        #   a_n = (1/N) sum_m b_m exp(-j 2 pi x_n cos(theta_m))
        #       = (1/N) sum_m b_m exp(-j 2 pi m n / N).
        n_idx = np.arange(N) - (N - 1) / 2.0
        amps = np.zeros(N, dtype=complex)
        for mi, bm in zip(m, b_m):
            if bm == 0.0:
                continue
            amps += bm * np.exp(-1j * 2.0 * pi * mi * n_idx / N) / N
        # for a target symmetric about broadside the excitations are real
        # (but can be negative, i.e. 180 deg phase)
        return amps

    def directivity(self, N, d_over_lambda, amps):
        # Exact pattern-integral directivity of the complex excitation.
        amps = np.asarray(amps)
        if not np.any(amps):
            return 0.0
        return self.directivity_integral(np.abs(amps), phases=np.angle(amps),
                                         d_over_lambda=d_over_lambda)

    # --- driver --------------------------------------------------------
    def woodward_lawson_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing
        shape = getattr(self.args, "shape", "flat_top")

        amps_c = self.amplitudes(N, d_over_lambda)
        mag = np.abs(amps_c)
        phase = np.angle(amps_c)                 # radians
        amps_norm = self.normalize(mag)
        # Report phases in (-180, 180] with numerical noise cleaned up, so a
        # real excitation prints as exactly 0 or 180 deg.
        phase_deg = np.round(np.degrees(phase), 9) + 0.0
        phase_deg[phase_deg <= -180.0] += 360.0
        phase_deg[mag < 1e-12 * max(mag.max(), 1e-300)] = 0.0
        positions = np.arange(N) * d_over_lambda
        D0 = self.directivity(N, d_over_lambda, amps_c)
        D0_db = 10.0 * log10(D0) if D0 > 0 else float("nan")

        if not self.args.variable_return:
            self.value_print("N", N)
            print("[*] Target shape =", shape)
            self.list_print("Amplitudes (normalized)", amps_norm)
            self.list_print("Excitation phase (deg)", phase_deg)
            self.report_spacing(d_over_lambda)
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        # The excitation phases already point the shaped beam at --scan, so
        # they are passed explicitly (no extra progressive steering phase).
        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, mag, positions=positions,
                                    phases=phase)
        if getattr(self.args, "plot", False):
            self.plot_pattern(mag, title="Woodward-Lawson (%s, N=%d)" % (shape, N),
                              style=getattr(self.args, "plot_style", "both"),
                              positions=positions, phases=phase)

        if self.args.variable_return:
            return amps_norm, phase_deg, D0_db
