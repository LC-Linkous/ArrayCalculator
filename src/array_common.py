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
from math import log10, cos, sin, sqrt, pi, radians, degrees
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

    # --- window sampling (Hann / Hamming / Blackman) -------------------
    def window_x(self, N, mode=None):
        # Normalized element coordinate x_n fed to the raised-cosine tapers,
        # with n = -(N-1)/2 .. (N-1)/2 centered on the array:
        #   'aperture'  (default): x = n / N. Samples a continuous window
        #       spanning the full aperture L = N d at the element centers, so
        #       the edge elements sit d/2 inside the window's end and keep a
        #       small nonzero excitation. This is the usual antenna convention
        #       (and makes Hann identical to cosine-squared).
        #   'symmetric': x = n / (N - 1). The DSP "symmetric" window: the edge
        #       elements land on the window's endpoints, so Hann (and, nearly,
        #       Blackman) drive them at zero -- effectively an (N-2)-element
        #       array with a wider beam and lower sidelobes.
        if mode is None:
            mode = getattr(self.args, "window_def", "aperture") or "aperture"
        n = np.arange(N) - (N - 1) / 2.0
        if mode == "symmetric":
            return n / (N - 1) if N > 1 else np.zeros(N)
        if mode != "aperture":
            raise ValueError("window_def must be 'aperture' or 'symmetric'")
        return n / N

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
        # Complex amplitudes are accepted: their phase adds to phi_n.
        amps = np.asarray(amps)
        N = len(amps)
        if d_over_lambda is None:
            d_over_lambda = getattr(self.args, "spacing", 0.5)
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
        # Complex weight per element, then the spatial phase per (angle,
        # element). Angles are processed in chunks so memory stays bounded
        # (~16 MB per chunk) for large N and fine sweeps.
        w = amps * np.exp(1j * phases)
        cos_t = np.cos(theta)
        chunk = max(1, (1 << 20) // max(N, 1))
        af = np.empty(cos_t.shape, dtype=float)
        for i in range(0, len(cos_t), chunk):
            arg = 2.0 * pi * positions[None, :] * cos_t[i:i + chunk, None]
            af[i:i + chunk] = np.abs(np.exp(1j * arg) @ w)
        return af

    def pattern_sweep(self, amps, n_points=721, positions=None, phases=None,
                      theta_range=(0.0, 180.0)):
        theta = np.linspace(theta_range[0], theta_range[1], n_points)
        af = self.array_factor(amps, theta, positions=positions, phases=phases)
        af_norm = af / af.max()
        af_db = 20.0 * np.log10(np.clip(af_norm, 1e-12, None))
        return theta, af_norm, af_db

    def resolve_geometry(self, N, positions=None, phases=None):
        # The (positions, phases) array_factor would use for these inputs:
        # equally spaced at args.spacing with the progressive --scan phase
        # unless given explicitly.
        if positions is None:
            positions = np.arange(N) * float(getattr(self.args, "spacing", 0.5))
        positions = np.asarray(positions, dtype=float)
        if phases is None:
            scan = getattr(self.args, "scan", 90.0)
            phases = -2.0 * pi * positions * cos(radians(scan))
        return positions, np.asarray(phases, dtype=float)

    # --- figures of merit (shared by every array) -----------------------
    def directivity_integral(self, amps, positions=None, phases=None,
                             d_over_lambda=None):
        # Exact array-factor directivity (isotropic elements) for any
        # positions, amplitudes and phases, including steered beams and
        # spacings with grating lobes:
        #   D = 2 |AF|max^2 / int_0^pi |AF|^2 sin(theta) dtheta
        # The denominator has a closed form. With w_n = a_n exp(j phi_n),
        #   int_0^pi |AF|^2 sin = sum_mn w_m w_n* 2 sinc(2 (x_m - x_n))
        # (numpy sinc), so D = |AF|max^2 / sum_mn w_m w_n* sinc(2 dx).
        amps = np.asarray(amps)
        if d_over_lambda is not None and positions is None:
            positions = np.arange(len(amps)) * float(d_over_lambda)
        positions, phases = self.resolve_geometry(len(amps), positions, phases)
        w = amps * np.exp(1j * phases)
        dx = positions[:, None] - positions[None, :]
        denom = np.real(np.sum(np.outer(w, np.conj(w)) * np.sinc(2.0 * dx)))
        if denom <= 0:
            return 0.0
        theta = np.linspace(0.0, 180.0, 40001)
        peak = np.max(self.array_factor(amps, theta, positions=positions,
                                        phases=phases)) ** 2
        return float(peak / denom)

    def hpbw_deg(self, amps, positions=None, phases=None, n_points=40001):
        # -3 dB width of the main beam, read from the pattern. AF depends on
        # theta only through cos(theta), so it continues past 0 and 180 deg
        # as a mirror image. Sweeping 180 deg either side of the peak means
        # the half-power walk never runs off the end, so endfire and
        # near-endfire beams get their full width in the plane cut. Returns
        # None if there is no half-power crossing (e.g. a single element).
        _, af0, _ = self.pattern_sweep(amps, n_points=n_points,
                                       positions=positions, phases=phases)
        peak_theta = 180.0 * np.argmax(af0) / (n_points - 1)
        theta, af, _ = self.pattern_sweep(
            amps, n_points=2 * n_points - 1, positions=positions, phases=phases,
            theta_range=(peak_theta - 180.0, peak_theta + 180.0))
        half = 1.0 / sqrt(2.0)
        peak = n_points - 1
        left = peak
        while left > 0 and af[left] > half:
            left -= 1
        right = peak
        while right < len(af) - 1 and af[right] > half:
            right += 1
        if af[left] > half or af[right] > half:
            return None

        # Linear interpolation between the samples bracketing each crossing.
        def cross(i_out, i_in):
            a0, a1 = af[i_out], af[i_in]
            t = (half - a0) / (a1 - a0) if a1 != a0 else 0.0
            return theta[i_out] + t * (theta[i_in] - theta[i_out])
        return float(cross(right, right - 1) - cross(left, left + 1))

    def hpbw_rad(self, amps, positions=None, phases=None):
        # hpbw_deg in radians (None if undefined), for the synthesis drivers.
        h = self.hpbw_deg(amps, positions=positions, phases=phases)
        return None if h is None else radians(h)

    def print_hpbw(self, hpbw_rad):
        if hpbw_rad is None:
            print("[*] HPBW = n/a (no half-power crossing)")
        else:
            self.value_print("HPBW", degrees(hpbw_rad), "deg")

    @staticmethod
    def hpbw_out(hpbw_rad):
        # Value returned by --variable_return: degrees, or None if undefined.
        return None if hpbw_rad is None else degrees(hpbw_rad)

    def closed_form_valid(self):
        # The textbook closed-form HPBW / directivity expressions used by the
        # binomial and Dolph-Tschebyscheff classes assume d = lambda/2 and a
        # broadside beam. Outside that, use the pattern-derived values.
        return (abs(float(getattr(self.args, "spacing", 0.5)) - 0.5) < 1e-9
                and abs(float(getattr(self.args, "scan", 90.0)) - 90.0) < 1e-9)

    @staticmethod
    def peak_sidelobe(theta, af_db, peak=None):
        # Highest lobe outside the main beam, in dB. The main beam is the
        # region between the nulls either side of the peak. Every other local
        # maximum counts, including grating lobes (which reach 0 dB) and lobes
        # sitting on the 0/180 deg ends of the sweep (the pattern is mirror-
        # symmetric there, so an end sample at least as high as its neighbour
        # is a lobe maximum). Returns -inf if there is no other lobe.
        af_db = np.asarray(af_db, dtype=float)
        n = len(af_db)
        if peak is None:
            peak = int(np.argmax(af_db))
        li = peak
        while li > 0 and af_db[li - 1] <= af_db[li]:
            li -= 1
        ri = peak
        while ri < n - 1 and af_db[ri + 1] <= af_db[ri]:
            ri += 1
        lobes = []
        for i in range(n):
            if li <= i <= ri:
                continue
            left_ok = i == 0 or af_db[i] >= af_db[i - 1]
            right_ok = i == n - 1 or af_db[i] >= af_db[i + 1]
            if left_ok and right_ok:
                lobes.append(af_db[i])
        return float(max(lobes)) if lobes else float("-inf")

    def warn_grating_lobes(self, amps, positions=None, phases=None):
        # Print a warning when a second lobe reaches within 3 dB of the main
        # beam (a grating lobe), since every figure of merit is then suspect.
        if getattr(self.args, "variable_return", False):
            return
        theta, _, af_db = self.pattern_sweep(amps, n_points=20001,
                                             positions=positions, phases=phases)
        if self.peak_sidelobe(theta, af_db) > -3.0:
            print("[!] Warning: grating lobe present (a second lobe within "
                  "3 dB of the main beam); reduce the spacing or the scan angle.")

    def export_pattern_csv(self, filename, amps, positions=None, phases=None,
                           n_points=721):
        theta, af_norm, af_db = self.pattern_sweep(amps, n_points=n_points,
                                                   positions=positions,
                                                   phases=phases)
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
        # theta is measured from the array axis, and a linear array's pattern
        # depends on it only through cos(theta), so the cut through any plane
        # containing the axis is mirror-symmetric about that axis. Draw the
        # swept 0..180 deg half on the top of the dial and its mirror image
        # (theta -> -theta) on the bottom: the array axis runs left-right
        # (0 deg on the right, 180 deg on the left) and broadside points up
        # and down. Angle labels show theta from the axis on both halves.
        th = np.radians(theta)
        ax.plot(th, af_db, color=self._ACCENT, lw=1.8, zorder=3)
        ax.plot(-th, af_db, color=self._ACCENT, lw=1.8, zorder=3)
        ax.set_theta_zero_location("E")
        ax.set_theta_direction(1)
        ticks = np.arange(0, 360, 30)
        ax.set_thetagrids(ticks, ["%d\u00b0" % (t if t <= 180 else 360 - t)
                                  for t in ticks])
        # dB labels along a diagonal, clear of the 0/180 deg angle labels and
        # of broadside/endfire beams. The rim is 0 dB, the center is `floor`;
        # neither needs its own label.
        ax.set_rlabel_position(22.5)
        ax.set_ylim(floor, 0)
        ax.set_rticks(list(range(floor, 1, 10)))
        ax.set_yticklabels(["" if r in (floor, 0) else str(r)
                            for r in range(floor, 1, 10)], fontsize=7)
        ax.grid(True, color=self._GRID, lw=0.6)
        if sll_db is not None:
            ring = np.linspace(0.0, 2.0 * pi, 361)
            ax.plot(ring, np.full_like(ring, -abs(sll_db)), color=self._SLLCOL,
                    ls="--", lw=1.0, label="SLL %g dB" % -abs(sll_db))
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

    def _element_axes(self, ax, positions, amps, phases=None):
        # Element layout along the array axis: stem height = |amplitude|,
        # x = element position in wavelengths. If any element carries a
        # nonzero excitation phase, color the markers by phase on a cyclic
        # colormap, with phases wrapped into [-180, 180) deg so a steered
        # array's progressive phase cycles through the colormap instead of
        # saturating at its ends. Dull for a uniform array; informative for
        # tapered, steered, or non-uniform geometries. Complex amplitudes are
        # split into magnitude (stem height) and phase (added to phases).
        amps = np.asarray(amps)
        if np.iscomplexobj(amps):
            extra = np.angle(amps)
            phases = extra if phases is None else np.asarray(phases, dtype=float) + extra
            amps = np.abs(amps)
        amps = np.asarray(amps, dtype=float)
        positions = np.asarray(positions, dtype=float)
        a = amps / amps.max() if amps.max() > 0 else amps
        ax.vlines(positions, 0, a, color=self._GRID, lw=1.2, zorder=1)
        ph = None
        if phases is not None:
            ph = np.degrees(np.asarray(phases, dtype=float))
            ph = (ph + 180.0) % 360.0 - 180.0          # wrap to [-180, 180)
        if ph is not None and np.max(np.abs(ph)) > 1e-6:
            sc = ax.scatter(positions, a, c=ph, cmap="twilight", s=38, zorder=3,
                            vmin=-180, vmax=180, edgecolor="k", linewidth=0.4)
            cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label("phase (deg)", fontsize=8)
            cb.ax.tick_params(labelsize=7)
        else:
            ax.scatter(positions, a, color=self._ACCENT, s=34, zorder=3,
                       edgecolor="k", linewidth=0.4)
        ax.set_xlabel("element position (wavelengths)")
        ax.set_ylabel("|amplitude|")
        ax.set_ylim(0, 1.08)
        if len(positions) > 1:
            pad = 0.3
            ax.set_xlim(positions.min() - pad, positions.max() + pad)
        ax.grid(True, color=self._GRID, lw=0.5, axis="y")
        ax.set_axisbelow(True)

    def plot_pattern(self, amps, title="Array Factor", sll_db=None,
                     style=None, floor=-60, savepath=None,
                     positions=None, phases=None, show_elements=None):
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
        if show_elements is None:
            show_elements = getattr(self.args, "elements_plot", False)
        # The plotted pattern must come from the same geometry the caller
        # describes: explicit positions/phases (evaluator, Woodward-Lawson)
        # when given, else equally spaced at --spacing with the --scan phase.
        theta, af_norm, af_db = self.pattern_sweep(amps, positions=positions,
                                                   phases=phases)
        af_db = np.clip(af_db, floor, 0)

        # Element panel uses the same resolved geometry.
        if show_elements:
            positions, phases = self.resolve_geometry(len(np.asarray(amps)),
                                                      positions, phases)

        # Figure layout: the pattern panel(s) on top, optional element row below.
        pattern_cols = 2 if style == "both" else 1
        if style == "polar":
            base = (5.4, 5.4)
        elif style == "rect":
            base = (6.6, 4.2)
        else:
            base = (10.6, 4.8)

        if show_elements:
            fig = plt.figure(figsize=(base[0], base[1] + 2.6))
            gs = fig.add_gridspec(2, pattern_cols, height_ratios=[base[1], 2.4])
            if style == "polar":
                self._polar_axes(fig.add_subplot(gs[0, 0], projection="polar"),
                                 theta, af_db, sll_db, floor)
            elif style == "rect":
                self._rect_axes(fig.add_subplot(gs[0, 0]), theta, af_db, sll_db, floor)
            else:
                self._polar_axes(fig.add_subplot(gs[0, 0], projection="polar"),
                                 theta, af_db, sll_db, floor)
                self._rect_axes(fig.add_subplot(gs[0, 1]), theta, af_db, sll_db, floor)
            # element panel spans the full width on the bottom row
            self._element_axes(fig.add_subplot(gs[1, :]), positions, amps, phases)
        elif style == "polar":
            fig = plt.figure(figsize=base)
            self._polar_axes(fig.add_subplot(111, projection="polar"),
                             theta, af_db, sll_db, floor)
        elif style == "rect":
            fig = plt.figure(figsize=base)
            self._rect_axes(fig.add_subplot(111), theta, af_db, sll_db, floor)
        else:  # both
            fig = plt.figure(figsize=base)
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