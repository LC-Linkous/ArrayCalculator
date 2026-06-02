#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator    array_common.py
#
#   Shared helpers for the linear-array synthesis classes
#   (BinomialArray, DolphTschebyscheff). Keeps the printing style
#   consistent with the rest of the calculator and holds the common
#   array-factor / pattern-export machinery so the two methods only
#   have to supply excitation amplitudes.
##--------------------------------------------------------------------\

import csv as _csv
from math import log10, cos, sin, pi, radians
import numpy as np
from pint import UnitRegistry  # pip install pint
ureg = UnitRegistry()

C_LIGHT = 3e8  # speed of light, matching the literal used elsewhere in the project


class ArrayCommon:
    def __init__(self, args):
        self.args = args

    # --- printing helpers ---------------------------------------------
    def unit_print(self, name, value, unit=None):
        # Length quantities, formatted via pint (mirrors the other classes).
        if unit is not None:
            print("[*]", name, "= {:.2f}".format((value * ureg.meter).to(self.args.unit)))
        else:
            print("[*]", name, "= {:.2f}".format((value * ureg.meter).to_compact()))

    def value_print(self, name, value, unit=None):
        # Dimensionless / non-length quantities (deg, dB, ratios, counts).
        if unit is not None:
            print("[*]", name, "= {:.2f} {}".format(value, unit))
        else:
            if isinstance(value, (int,)) or float(value).is_integer():
                print("[*]", name, "=", int(value))
            else:
                print("[*]", name, "= {:.2f}".format(value))

    def list_print(self, name, values):
        formatted = ", ".join("{:.3f}".format(v) for v in values)
        print("[*]", name, "= [" + formatted + "]")

    # --- frequency / spacing ------------------------------------------
    def wavelength(self):
        if getattr(self.args, "frequency", None):
            return C_LIGHT / self.args.frequency
        return None

    def report_spacing(self, d_over_lambda):
        # If a frequency was supplied, print wavelength and physical spacing.
        lam = self.wavelength()
        if lam is not None:
            self.unit_print("Wavelength", lam, self.args.unit)
            self.unit_print("Element spacing d", d_over_lambda * lam, self.args.unit)
        else:
            self.value_print("Element spacing d", d_over_lambda, "lambda")

    # --- array factor / pattern ---------------------------------------
    def array_factor(self, amps, theta_deg, d_over_lambda=None, scan_deg=None):
        # |AF(theta)| for a linear array of equally spaced elements.
        # psi = 2*pi*(d/lambda)*cos(theta) - beta, with beta the progressive
        # phase that steers the beam to scan_deg (broadside => scan_deg=90).
        amps = np.asarray(amps, dtype=float)
        N = len(amps)
        if d_over_lambda is None:
            d_over_lambda = self.args.spacing
        if scan_deg is None:
            scan_deg = getattr(self.args, "scan", 90.0)

        n = np.arange(N)
        beta = -2.0 * pi * d_over_lambda * cos(radians(scan_deg))
        theta = np.radians(np.asarray(theta_deg, dtype=float))
        psi = 2.0 * pi * d_over_lambda * np.cos(theta)[:, None] + beta
        af = np.abs(np.sum(amps[None, :] * np.exp(1j * n[None, :] * psi), axis=1))
        return af

    def pattern_sweep(self, amps, n_points=721):
        theta = np.linspace(0.0, 180.0, n_points)
        af = self.array_factor(amps, theta)
        af_norm = af / af.max()
        af_db = 20.0 * np.log10(np.clip(af_norm, 1e-12, None))
        return theta, af_norm, af_db

    def export_pattern_csv(self, filename, amps):
        theta, af_norm, af_db = self.pattern_sweep(amps)
        with open(filename, "w", newline="") as f:
            w = _csv.writer(f)
            w.writerow(["theta_deg", "AF_linear", "AF_dB"])
            for t, a, d in zip(theta, af_norm, af_db):
                w.writerow(["{:.3f}".format(t), "{:.6f}".format(a), "{:.4f}".format(d)])
        print("[*] Pattern CSV generated: " + filename)

    # --- plotting (future expansion; matplotlib imported lazily) -------
    def plot_pattern(self, amps, title="Array Factor", sll_db=None):
        # Lazy import keeps matplotlib an optional dependency: the core
        # calculator runs without it installed unless --plot is used.
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("[*] matplotlib not installed; skipping plot. "
                  "Install it or use --csv for pattern data.")
            return
        theta, af_norm, af_db = self.pattern_sweep(amps)
        af_db = np.clip(af_db, -60, 0)
        ax = plt.subplot(111, projection="polar")
        ax.plot(np.radians(theta), af_db)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rlabel_position(90)
        ax.set_ylim(-60, 0)
        if sll_db is not None:
            ax.plot(np.radians(theta), np.full_like(theta, -abs(sll_db)),
                    linestyle="--", linewidth=0.8)
        ax.set_title(title)
        plt.show()