#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  array_evaluator.py
#
#   Evaluator (not a synthesis method): given an arbitrary array geometry
#   -- element POSITIONS (in wavelengths), and optionally per-element
#   AMPLITUDES and PHASES -- it runs the full analysis pipeline (pattern
#   sweep, HPBW, peak sidelobe level, directivity-by-integral) and can
#   export/plot the result, exactly like the synthesis classes.
#
#   This is the boundary between this calculator and an external optimizer:
#   the optimizer PROPOSES positions/amplitudes/phases, the evaluator SCORES
#   them. No optimizer lives here. Two entry points:
#
#     * evaluate(positions, amplitudes=None, phases=None) -> dict
#         A pure method for programmatic use (the optimizer's fitness call).
#         amplitudes default to uniform, phases to zero (broadside, no steer).
#
#     * CLI subcommand `evaluate` that reads a geometry CSV with columns
#         position_lambda[, amplitude][, phase_deg]
#       and prints / plots / re-exports the analysis.
#
#   Directivity is the array-factor (isotropic-element) value from the
#   pattern integral, valid for arbitrary (non-uniform) geometries; the
#   simple 2 N d / lambda closed form does not apply once spacing is
#   unequal, so it is not used here.
##--------------------------------------------------------------------\

import csv as _csv
from math import log10, sqrt, pi, radians, degrees
import numpy as np

from array_common import ArrayCommon


class ArrayEvaluator(ArrayCommon):
    def __init__(self, args=None):
        super().__init__(args)
        self.args = args

    # --- core scoring --------------------------------------------------
    def evaluate(self, positions, amplitudes=None, phases=None,
                 n_points=40001):
        """Score an arbitrary linear-array geometry.

        positions : element positions in wavelengths (array-like, length N)
        amplitudes: per-element excitation magnitude (default all ones)
        phases    : per-element excitation phase in RADIANS (default zeros)

        Returns a dict with the realized pattern and its key figures:
          N, positions, amplitudes, phases,
          theta_deg, af_norm, af_db,
          peak_theta_deg, hpbw_deg (None if no clean single peak),
          peak_sidelobe_db, directivity, directivity_db
        """
        positions = np.asarray(positions, dtype=float)
        N = len(positions)
        if amplitudes is None:
            amplitudes = np.ones(N)
        amplitudes = np.asarray(amplitudes, dtype=float)
        if phases is None:
            phases = np.zeros(N)
        phases = np.asarray(phases, dtype=float)

        theta = np.linspace(0.0, 180.0, n_points)
        af = self.array_factor(amplitudes, theta, positions=positions,
                               phases=phases)
        af_norm = af / af.max()
        af_db = 20.0 * np.log10(np.clip(af_norm, 1e-12, None))

        peak = int(np.argmax(af_norm))
        peak_theta = float(theta[peak])
        hpbw = self._hpbw_from_pattern(theta, af_norm, peak)
        sll = self._peak_sidelobe(theta, af_db, peak)
        D0 = self._directivity_integral(amplitudes, positions, phases)

        return {
            "N": N,
            "positions": positions,
            "amplitudes": amplitudes,
            "phases": phases,
            "theta_deg": theta,
            "af_norm": af_norm,
            "af_db": af_db,
            "peak_theta_deg": peak_theta,
            "hpbw_deg": (degrees(hpbw) if hpbw is not None else None),
            "peak_sidelobe_db": sll,
            "directivity": D0,
            "directivity_db": (10.0 * log10(D0) if D0 > 0 else float("nan")),
        }

    # --- helpers -------------------------------------------------------
    def _hpbw_from_pattern(self, theta, af_norm, peak):
        half = 1.0 / sqrt(2.0)
        left = peak
        while left > 0 and af_norm[left] > half:
            left -= 1
        right = peak
        while right < len(af_norm) - 1 and af_norm[right] > half:
            right += 1
        if left == 0 or right == len(af_norm) - 1:
            return None  # beam runs off the visible region; HPBW ill-defined
        return radians(theta[right] - theta[left])

    def _peak_sidelobe(self, theta, af_db, peak):
        # nearest null on each side of the main beam, then the highest local
        # maximum outside that main-beam region.
        n = len(af_db)
        li = peak
        while li > 0 and af_db[li - 1] <= af_db[li]:
            li -= 1
        ri = peak
        while ri < n - 1 and af_db[ri + 1] <= af_db[ri]:
            ri += 1
        peaks = [af_db[i] for i in range(1, n - 1)
                 if (i < li or i > ri)
                 and af_db[i] >= af_db[i - 1] and af_db[i] >= af_db[i + 1]
                 and af_db[i] < -0.5]
        return float(max(peaks)) if peaks else float("-inf")

    def _directivity_integral(self, amps, positions, phases, n=200001):
        # D = 2 |AF(peak)|^2 / integral_0^pi |AF|^2 sin(theta) dtheta,
        # isotropic elements. Valid for arbitrary positions/phases.
        th = np.linspace(1e-6, np.pi - 1e-6, n)
        arg = (2.0 * pi * positions[None, :] * np.cos(th)[:, None]
               + phases[None, :])
        af = np.abs((amps[None, :] * np.exp(1j * arg)).sum(axis=1))
        p = af ** 2
        peak_power = p.max()
        integral = np.trapezoid(p * np.sin(th), th)
        return 2.0 * peak_power / integral if integral > 0 else 0.0

    # --- CSV geometry I/O ---------------------------------------------
    @staticmethod
    def read_geometry_csv(filename):
        """Read a geometry CSV: position_lambda[, amplitude][, phase_deg].

        amplitude defaults to 1.0, phase_deg to 0.0 when columns are absent.
        Returns (positions, amplitudes, phases_radians).
        """
        positions, amps, phases = [], [], []
        with open(filename, newline="") as f:
            reader = _csv.DictReader(f)
            cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
            pcol = cols.get("position_lambda") or cols.get("position")
            acol = cols.get("amplitude")
            qcol = cols.get("phase_deg") or cols.get("phase")
            if pcol is None:
                raise ValueError("geometry CSV needs a 'position_lambda' column")
            for row in reader:
                positions.append(float(row[pcol]))
                amps.append(float(row[acol]) if acol and row.get(acol) not in (None, "") else 1.0)
                phases.append(radians(float(row[qcol])) if qcol and row.get(qcol) not in (None, "") else 0.0)
        return (np.array(positions), np.array(amps), np.array(phases))

    # --- driver (CLI) --------------------------------------------------
    def evaluate_calculator(self):
        geom = self.args.geometry
        positions, amps, phases = self.read_geometry_csv(geom)
        res = self.evaluate(positions, amps, phases)

        if not self.args.variable_return:
            self.value_print("N", res["N"])
            self.list_print("Positions (lambda)", positions)
            self.list_print("Amplitudes", amps)
            if np.any(phases):
                self.list_print("Phases (deg)", np.degrees(phases))
            self.value_print("Main beam at", res["peak_theta_deg"], "deg")
            if res["hpbw_deg"] is not None:
                self.value_print("HPBW", res["hpbw_deg"], "deg")
            else:
                print("[*] HPBW = n/a (beam not a single clean peak in [0,180])")
            self.value_print("Peak sidelobe", res["peak_sidelobe_db"], "dB")
            self.value_print("Directivity", res["directivity"])
            self.value_print("Directivity", res["directivity_db"], "dB")

        if getattr(self.args, "csv", None):
            theta, af_norm, af_db = (res["theta_deg"], res["af_norm"], res["af_db"])
            with open(self.args.csv, "w", newline="") as f:
                w = _csv.writer(f)
                w.writerow(["theta_deg", "AF_linear", "AF_dB"])
                for t, a, d in zip(theta, af_norm, af_db):
                    w.writerow(["{:.3f}".format(t), "{:.6f}".format(a), "{:.4f}".format(d)])
            print("[*] Pattern CSV generated: " + self.args.csv)
        if getattr(self.args, "plot", False):
            self.plot_pattern(amps, title="Evaluated geometry (N=%d)" % res["N"],
                              style=getattr(self.args, "plot_style", "both"))

        if self.args.variable_return:
            return res