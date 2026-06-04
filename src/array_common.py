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

    # --- normalization ------------------------------------------------
    def normalize(self, amps, mode=None, eps=1e-9):
        # Scale amplitudes by a single positive factor; never reorders or
        # reshapes. 'center' -> peak element becomes 1.0; 'edge' -> smallest
        # significant element becomes 1.0. Edge mode guards against the
        # near-zero edge samples some tapers produce (e.g. Blackman), which
        # would otherwise divide by ~0. Both modes describe the same physical
        # array and produce an identical radiation pattern; only the printed
        # scale differs.
        amps = np.asarray(amps, dtype=float)
        if mode is None:
            mode = getattr(self.args, "norm", "center")
        if mode == "edge":
            sig = amps[amps > eps]
            ref = sig.min() if sig.size else amps.max()
        else:  # center
            ref = amps.max()
        return amps / ref

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
    def array_factor(self, amps, theta_deg, d_over_lambda=None, scan_deg=None,
                     positions=None, phases=None):
        # |AF(theta)| for a linear array along one axis. The general form is
        #   AF(theta) = sum_n a_n exp( j [ 2 pi x_n cos(theta) + phi_n ] )
        # where x_n is the element position in wavelengths and phi_n is the
        # per-element excitation phase (radians).
        #
        # Backward-compatible defaults reproduce the equally-spaced, linearly
        # -steered case exactly: if positions is None the elements sit at
        # x_n = n * (d/lambda); if phases is None a progressive steering phase
        # phi_n = -2 pi x_n cos(scan) points the beam at scan_deg (broadside =>
        # scan_deg = 90 => phi_n = 0). Pass positions and/or phases explicitly
        # to evaluate an arbitrary (non-uniform, individually phased) geometry.
        amps = np.asarray(amps, dtype=float)
        N = len(amps)
        if d_over_lambda is None:
            d_over_lambda = self.args.spacing
        if scan_deg is None:
            scan_deg = getattr(self.args, "scan", 90.0)

        if positions is None:
            positions = np.arange(N) * float(d_over_lambda)
        else:
            positions = np.asarray(positions, dtype=float)
        if phases is None:
            # progressive phase that steers a uniform array to scan_deg
            phases = -2.0 * pi * positions * cos(radians(scan_deg))
        else:
            phases = np.asarray(phases, dtype=float)

        theta = np.radians(np.asarray(theta_deg, dtype=float))
        # spatial phase per (angle, element) plus the per-element excitation phase
        arg = (2.0 * pi * positions[None, :] * np.cos(theta)[:, None]
               + phases[None, :])
        af = np.abs(np.sum(amps[None, :] * np.exp(1j * arg), axis=1))
        return af

    def pattern_sweep(self, amps, n_points=721, positions=None, phases=None):
        theta = np.linspace(0.0, 180.0, n_points)
        af = self.array_factor(amps, theta, positions=positions, phases=phases)
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

    # --- plotting (matplotlib imported lazily) ------------------------
    # Color/style constants shared by both plot styles.
    _ACCENT = "#1f6feb"   # main beam trace
    _SLLCOL = "#d1495b"   # sidelobe-level reference
    _GRID = "#c9ced6"     # grid lines

    def _polar_axes(self, ax, theta, af_db, sll_db, floor):
        # A linear array's pattern is swept over theta in [0, 180] and is
        # symmetric about the boresight axis. We mirror the sweep to the
        # lower half so the dial shows the true symmetry with the beam
        # pointing up (90 deg physical -> top), and only draw the
        # meaningful 0..180 half-plane rather than a full empty circle.
        th = np.radians(theta)
        ax.plot(th, af_db, color=self._ACCENT, lw=1.8, zorder=3)
        ax.plot(-th, af_db, color=self._ACCENT, lw=1.8, zorder=3)
        ax.set_theta_zero_location("E")
        ax.set_theta_direction(1)
        ax.set_thetamin(0)
        ax.set_thetamax(180)
        ax.set_rlabel_position(0)          # dB labels along the 0 deg axis
        ax.set_ylim(floor, 0)
        ax.set_rticks(list(range(floor, 1, 10)))
        ax.grid(True, color=self._GRID, lw=0.6)
        if sll_db is not None:
            lbl = "SLL %g dB" % -abs(sll_db)
            ax.plot(th, np.full_like(th, -abs(sll_db)), color=self._SLLCOL,
                    ls="--", lw=1.0, label=lbl)
            ax.plot(-th, np.full_like(th, -abs(sll_db)),
                    color=self._SLLCOL, ls="--", lw=1.0)
            ax.legend(loc="lower center", fontsize=8,
                      bbox_to_anchor=(0.5, -0.12), frameon=False)

    def _rect_axes(self, ax, theta, af_db, sll_db, floor):
        ax.plot(theta, af_db, color=self._ACCENT, lw=1.8, zorder=3, label="AF")
        ax.set_xlim(0, 180)
        ax.set_ylim(floor, 1)
        ax.set_xticks(range(0, 181, 30))
        ax.set_xlabel("theta (deg)")
        ax.set_ylabel("normalized |AF| (dB)")
        ax.grid(True, color=self._GRID, lw=0.6)
        ax.axhline(0, color="#888888", lw=0.6)
        if sll_db is not None:
            ax.axhline(-abs(sll_db), color=self._SLLCOL, ls="--", lw=1.0,
                       label="SLL %g dB" % -abs(sll_db))
        ax.legend(loc="upper right", fontsize=8)

    def plot_pattern(self, amps, title="Array Factor", sll_db=None,
                     style=None, floor=-60, savepath=None):
        # Lazy import keeps matplotlib an optional dependency: the core
        # calculator runs without it installed unless --plot is used.
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("[*] matplotlib not installed; skipping plot. "
                  "Install it or use --csv for pattern data.")
            return
        # Fall back to the CLI flags when the caller doesn't override them,
        # so --plot-style and --save work for every array without each
        # driver having to forward them explicitly.
        if style is None:
            style = getattr(self.args, "plot_style", "both") or "both"
        if savepath is None:
            savepath = getattr(self.args, "save", None)
        theta, af_norm, af_db = self.pattern_sweep(amps)
        af_db = np.clip(af_db, floor, 0)

        if style == "polar":
            fig = plt.figure(figsize=(5.4, 5.4))
            self._polar_axes(fig.add_subplot(111, projection="polar"),
                             theta, af_db, sll_db, floor)
        elif style == "rect":
            fig = plt.figure(figsize=(6.6, 4.2))
            self._rect_axes(fig.add_subplot(111), theta, af_db, sll_db, floor)
        else:  # both
            fig = plt.figure(figsize=(10.6, 4.8))
            self._polar_axes(fig.add_subplot(121, projection="polar"),
                             theta, af_db, sll_db, floor)
            self._rect_axes(fig.add_subplot(122), theta, af_db, sll_db, floor)

        fig.suptitle(title, fontsize=12)
        fig.tight_layout()
        if savepath is not None:
            fig.savefig(savepath, dpi=110, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()