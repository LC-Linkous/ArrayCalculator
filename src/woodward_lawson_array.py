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
#     * flat-top  : uniform amplitude across a sector, zero outside
#     * cosecant^2: ~csc^2(theta) over a fill sector (classic radar
#                   ground-coverage pattern), with a floor outside
#
#   Reference: Woodward & Lawson (1948); Balanis, Antenna Theory.
##--------------------------------------------------------------------\

from math import log10, sqrt, pi, radians, degrees
import numpy as np

from array_common import ArrayCommon


class WoodwardLawsonArray(ArrayCommon):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    # --- desired pattern targets --------------------------------------
    def _target(self, theta_deg):
        # Desired |pattern| as a function of theta (deg). Centered on the
        # scan angle; the fill sector half-width and floor are CLI options.
        shape = getattr(self.args, "shape", "flat_top")
        center = getattr(self.args, "scan", 90.0)
        half = getattr(self.args, "sector", 30.0)      # half-width in deg
        floor = getattr(self.args, "floor", 0.0)        # linear floor outside
        th = np.asarray(theta_deg, dtype=float)
        out = np.full_like(th, floor)
        lo, hi = center - half, center + half
        inside = (th >= lo) & (th <= hi)
        if shape == "flat_top":
            out[inside] = 1.0
        elif shape == "cosecant_squared":
            # csc^2 coverage: amplitude ~ csc(theta)/csc(theta_ref), clipped.
            # Defined relative to the near edge of the sector so it peaks at
            # the lower angle and falls across the fill region.
            ref = np.radians(np.clip(lo, 1.0, 179.0))
            thr = np.radians(np.clip(th[inside], 1.0, 179.0))
            val = np.sin(ref) / np.sin(thr)            # csc(th)/csc(ref)
            out[inside] = np.clip(val, 0.0, 1.0)
        else:
            raise ValueError("unknown shape %r" % shape)
        return out

    # --- synthesis -----------------------------------------------------
    def amplitudes(self, N, d_over_lambda=None):
        # Composing-beam sample points: for an N-element array there are N
        # orthogonal beam positions at psi_m = 2 pi m / N, i.e. the sin-space
        # sample grid. We sample the desired pattern at the corresponding
        # angles and synthesize excitations as the inverse-DFT superposition.
        if d_over_lambda is None:
            d_over_lambda = self.args.spacing
        m = np.arange(N) - (N - 1) / 2.0
        # orthogonal sample directions in cos(theta) space:
        #   cos(theta_m) = m / (N * d/lambda)
        cos_tm = m / (N * d_over_lambda)
        valid = np.abs(cos_tm) <= 1.0                  # only real (visible) angles
        theta_m = np.degrees(np.arccos(np.clip(cos_tm, -1.0, 1.0)))
        b_m = np.where(valid, self._target(theta_m), 0.0)   # desired samples

        # Element excitations: a_n = sum_m b_m exp(-j 2 pi (m_idx) (n) / N).
        # Build directly from the composing-beam superposition.
        n_idx = np.arange(N) - (N - 1) / 2.0
        amps = np.zeros(N, dtype=complex)
        for mi, bm in zip(m, b_m):
            if bm == 0.0:
                continue
            amps += bm * np.exp(1j * 2.0 * pi * mi * n_idx / N) / N
        # for a symmetric real target the excitations are (near) real
        return amps

    def hpbw_or_none(self, amps):
        # Shaped beams do not have a single meaningful HPBW; report the -3 dB
        # width of the main feature only if one clean peak exists, else None.
        mag = np.abs(amps)
        theta, af_norm, _ = self.pattern_sweep(mag, n_points=40001)
        return None  # shaped beam: HPBW not reported

    def directivity(self, N, d_over_lambda, amps):
        mag = np.abs(np.asarray(amps))
        eff = (mag.sum() ** 2) / (len(mag) * np.sum(mag ** 2)) if np.any(mag) else 0.0
        return eff * 2.0 * N * d_over_lambda

    # --- driver --------------------------------------------------------
    def woodward_lawson_array_calculator(self):
        N = self.args.elements
        d_over_lambda = self.args.spacing
        shape = getattr(self.args, "shape", "flat_top")

        amps_c = self.amplitudes(N, d_over_lambda)
        mag = np.abs(amps_c)
        amps_norm = self.normalize(mag)
        phase_deg = np.degrees(np.angle(amps_c))
        D0 = self.directivity(N, d_over_lambda, mag)
        D0_db = 10.0 * log10(D0) if D0 > 0 else float("nan")

        if not self.args.variable_return:
            self.value_print("N", N)
            print("[*] Target shape =", shape)
            self.list_print("Amplitudes (normalized)", amps_norm)
            if self.args.verbose:
                self.list_print("Excitation phase (deg)", phase_deg)
            self.report_spacing(d_over_lambda)
            self.value_print("Directivity", D0)
            self.value_print("Directivity", D0_db, "dB")

        if getattr(self.args, "csv", None):
            self.export_pattern_csv(self.args.csv, mag)
        if getattr(self.args, "plot", False):
            self.plot_pattern(mag, title="Woodward-Lawson (%s, N=%d)" % (shape, N),
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return amps_norm, phase_deg, D0_db
