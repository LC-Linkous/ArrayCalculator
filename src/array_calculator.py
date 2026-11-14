#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  array_calculator.py
#
#   This is the main entry point of the array calculator.
#   This project is designed as an AntennaCAT compatible calcuator.
#   Refer to the README for argument formats.
##--------------------------------------------------------------------\


import argparse

__version__ = "1.0.0"   # keep in step with pyproject.toml

from binomial_array import BinomialArray
from dolph_tschebyscheff import DolphTschebyscheff
from uniform_array import UniformArray
from taylor_array import TaylorArray
from triangular_array import TriangularArray
from cosine_array import CosineArray
from cosine_squared_array import CosineSquaredArray
from hann_array import HannArray
from hamming_array import HammingArray
from blackman_array import BlackmanArray
from bartlett_array import BartlettArray
from kaiser_array import KaiserArray
from villeneuve_array import VilleneuveArray
from woodward_lawson_array import WoodwardLawsonArray
from array_evaluator import ArrayEvaluator


# Arrays whose only inputs are the common ones (N, spacing, freq, plotting,
# etc.) -- i.e. no method-specific flags like -sll. Maps the subcommand name
# to (class, driver-method-name). The parametric arrays (Dolph, Taylor,
# Villeneuve, Kaiser, Woodward-Lawson, and the evaluator) are built separately
# below because they add their own arguments.
SIMPLE_ARRAYS = {
    "uniform_array": (UniformArray, "uniform_array_calculator"),
    "binomial_array": (BinomialArray, "binomial_array_calculator"),
    "triangular_array": (TriangularArray, "triangular_array_calculator"),
    "bartlett_array": (BartlettArray, "bartlett_array_calculator"),
    "cosine_array": (CosineArray, "cosine_array_calculator"),
    "cosine_squared_array": (CosineSquaredArray, "cosine_squared_array_calculator"),
    "hann_array": (HannArray, "hann_array_calculator"),
    "hamming_array": (HammingArray, "hamming_array_calculator"),
    "blackman_array": (BlackmanArray, "blackman_array_calculator"),
}


def _add_common_args(p, synthesis=True):
    # The argument set every array subcommand shares. Centralized here so
    # all arrays stay in parity (spacing, steering, units, output, plotting)
    # and new flags only have to be added in one place. synthesis=False
    # (the evaluator) leaves out the geometry/normalization flags, because
    # the geometry CSV supplies positions and phases; accepting them and
    # silently ignoring them would mislead.
    p.add_argument('--help', action='help', default=argparse.SUPPRESS,
                   help='Show this help message and exit')
    p.add_argument('--verbose', action='store_true')
    if synthesis:
        _add_geometry_args(p)
    p.add_argument('--csv', type=str, required=False,
                   help='Export radiation pattern to CSV file')
    p.add_argument('--plot', action='store_true', required=False, default=False,
                   help='Plot radiation pattern (requires matplotlib)')
    p.add_argument('--plot-style', dest='plot_style', type=str,
                   choices=['polar', 'rect', 'both'], required=False, default='both',
                   help='Plot layout: polar dial, rectangular dB, or both (default both)')
    p.add_argument('--show-elements', dest='elements_plot', action='store_true',
                   required=False, default=False,
                   help='Add an element-layout panel (amplitude as stem height, '
                        'phase as marker color) beneath the pattern plot')
    p.add_argument('--save', type=str, required=False, default=None,
                   help='Save the plot to this path instead of opening a window '
                        '(implies --plot; works without a display)')
    p.add_argument('--variable_return', action='store_true', required=False, default=False,
                   help='Return Variables instead of printing')


def _add_geometry_args(p):
    # Flags that describe a synthesized, equally spaced array.
    p.add_argument('-N', '--elements', type=int, required=True,
                   help='Number of array elements (>= 2)')
    p.add_argument('-f', '--frequency', type=float, required=False,
                   help='Frequency in Hz (optional; gives physical spacing)')
    p.add_argument('-d', '--spacing', type=float, required=False, default=0.5,
                   help='Element spacing as a fraction of wavelength (default 0.5)')
    p.add_argument('--scan', type=float, required=False, default=90.0,
                   help='Beam steering angle in degrees (default 90 = broadside)')
    p.add_argument('-u', '--unit', type=str,
                   choices=['meter', 'centimeter', 'millimeter', 'inch'],
                   required=False, default='centimeter',
                   help='Unit of measurement')
    p.add_argument('--norm', type=str, choices=['edge', 'center'],
                   required=False, default='center',
                   help='Normalize amplitudes to the edge or center element '
                        '(default center; Dolph-Tschebyscheff conventionally uses edge)')


def _validate(parser, args):
    # Reject inputs the synthesis math cannot handle, with a clear message
    # instead of a NaN, ZeroDivisionError or silently corrupted result.
    name = args.subparser_name

    def need(ok, msg):
        if not ok:
            parser.error(msg)

    if name == 'evaluate':
        return
    N = args.elements
    need(N >= 2, "-N/--elements must be at least 2 (got %d)" % N)
    need(args.spacing > 0, "-d/--spacing must be > 0 (got %g)" % args.spacing)
    need(0.0 <= args.scan <= 180.0, "--scan must be within 0..180 deg (got %g)" % args.scan)
    if args.frequency is not None:
        need(args.frequency > 0, "-f/--frequency must be > 0 (got %g)" % args.frequency)
    if hasattr(args, 'sidelobe_level'):
        need(args.sidelobe_level > 0,
             "-sll must be a positive dB value, e.g. 26 (got %g)" % args.sidelobe_level)
    if hasattr(args, 'nbar'):
        need(args.nbar >= 1, "-nbar must be >= 1 (got %d)" % args.nbar)
    if name == 'villeneuve_array':
        nmax = (N + 1) // 2
        need(args.nbar <= nmax,
             "-nbar must be <= (N+1)//2 = %d for N=%d (got %d): the inner "
             "zeros would overlap their mirror images" % (nmax, N, args.nbar))
    if name == 'kaiser_array':
        need(args.beta >= 0, "-beta must be >= 0 (got %g)" % args.beta)
    if name == 'woodward_lawson_array':
        need(args.sector > 0, "--sector must be > 0 deg (got %g)" % args.sector)
        need(args.floor >= 0, "--floor must be >= 0 (got %g)" % args.floor)
        if args.shape == 'cosecant_squared':
            need(abs(args.scan - 90.0) > 1e-9,
                 "cosecant_squared needs --scan away from 90 (the beam peak "
                 "must sit off the horizon, e.g. --scan 80)")


class ArrayCalculator:
    def __init__(self, a=None):
        main_parser = argparse.ArgumentParser(description='Array Calculator', add_help=False)
        main_parser.add_argument('--help', action='help', default=argparse.SUPPRESS,
                                 help='Show this help message and exit')
        main_parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

        subparsers = main_parser.add_subparsers(help='sub-command help',
                                                dest='subparser_name')

        # Arrays that take only the common arguments.
        for name in SIMPLE_ARRAYS:
            sp = subparsers.add_parser(name, add_help=False)
            _add_common_args(sp)

        # Raised-cosine tapers also take --window-def (sampling convention).
        for name in ('hann_array', 'hamming_array', 'blackman_array'):
            subparsers.choices[name].add_argument(
                '--window-def', dest='window_def', type=str,
                choices=['aperture', 'symmetric'], default='aperture',
                help='Taper sampling: aperture = n/N (default; edge elements '
                     'stay excited), symmetric = n/(N-1) (DSP window; Hann '
                     'edge elements go to zero)')

        # DOLPH-TSCHEBYSCHEFF ARRAY (adds -sll and --norm)
        dt = subparsers.add_parser('dolph_tschebyscheff', add_help=False)
        _add_common_args(dt)
        dt.add_argument('-sll', '--sidelobe_level', type=float, required=True,
                        help='Desired sidelobe level in dB (e.g. 26, 30)')
        # Dolph conventionally reports edge-normalized coefficients; override
        # the common default (center) for this subcommand only.
        dt.set_defaults(norm='edge')

        # TAYLOR ARRAY (adds -sll and -nbar)
        ta = subparsers.add_parser('taylor_array', add_help=False)
        _add_common_args(ta)
        ta.add_argument('-sll', '--sidelobe_level', type=float, required=True,
                        help='Desired (near-in) sidelobe level in dB (e.g. 25, 30)')
        ta.add_argument('-nbar', '--nbar', type=int, required=False, default=5,
                        help='Number of near-in sidelobes held near the design level '
                             '(default 5)')

        # VILLENEUVE ARRAY (discrete Taylor; adds -sll and -nbar)
        vl = subparsers.add_parser('villeneuve_array', add_help=False)
        _add_common_args(vl)
        vl.add_argument('-sll', '--sidelobe_level', type=float, required=True,
                        help='Desired (near-in) sidelobe level in dB')
        vl.add_argument('-nbar', '--nbar', type=int, required=False, default=5,
                        help='Number of near-in sidelobes held near the design level '
                             '(default 5)')

        # KAISER ARRAY (parametric taper; adds -beta)
        ka = subparsers.add_parser('kaiser_array', add_help=False)
        _add_common_args(ka)
        ka.add_argument('-beta', '--beta', type=float, required=False, default=6.0,
                        help='Kaiser shape parameter (larger = lower sidelobes, '
                             'wider beam; default 6.0)')

        # WOODWARD-LAWSON ARRAY (shaped beam; adds --shape, --sector, --floor)
        wl = subparsers.add_parser('woodward_lawson_array', add_help=False)
        _add_common_args(wl)
        wl.add_argument('--shape', type=str, choices=['flat_top', 'cosecant_squared'],
                        required=False, default='flat_top',
                        help='Target pattern shape (default flat_top)')
        wl.add_argument('--sector', type=float, required=False, default=30.0,
                        help='Half-width of the shaped sector in degrees (default 30)')
        wl.add_argument('--floor', type=float, required=False, default=0.0,
                        help='Desired linear amplitude outside the sector (default 0)')

        # EVALUATOR (scores an arbitrary geometry from a CSV; not a synthesis)
        ev = subparsers.add_parser('evaluate', add_help=False)
        _add_common_args(ev, synthesis=False)
        ev.add_argument('-g', '--geometry', type=str, required=True,
                        help='Geometry CSV: position_lambda[, amplitude][, phase_deg]')

        self.args = main_parser.parse_args(a)
        if self.args.subparser_name is None:
            main_parser.error("choose an array subcommand (see --help)")
        _validate(subparsers.choices[self.args.subparser_name], self.args)
        # --save implies --plot so the user doesn't have to pass both.
        if getattr(self.args, 'save', None):
            self.args.plot = True
        self.calcedParams = None  # to catch returned vars if they exist

    def main(self, args):
        name = args.subparser_name

        if name in SIMPLE_ARRAYS:
            cls, driver = SIMPLE_ARRAYS[name]
            obj = cls(args)
            self.calcedParams = getattr(obj, driver)()

        elif name == 'dolph_tschebyscheff':
            dt = DolphTschebyscheff(args)
            self.calcedParams = dt.dolph_tschebyscheff_calculator()

        elif name == 'taylor_array':
            t = TaylorArray(args)
            self.calcedParams = t.taylor_array_calculator()

        elif name == 'villeneuve_array':
            v = VilleneuveArray(args)
            self.calcedParams = v.villeneuve_array_calculator()

        elif name == 'kaiser_array':
            k = KaiserArray(args)
            self.calcedParams = k.kaiser_array_calculator()

        elif name == 'woodward_lawson_array':
            w = WoodwardLawsonArray(args)
            self.calcedParams = w.woodward_lawson_array_calculator()

        elif name == 'evaluate':
            e = ArrayEvaluator(args)
            self.calcedParams = e.evaluate_calculator()

    def getArgs(self):
        return self.args

    def getCalcedParams(self):
        return self.calcedParams


def main(argv=None):
    # Console-script entry point (`arraycalculator ...`) and `python
    # array_calculator.py ...`. Returns the --variable_return result, if any.
    shell = ArrayCalculator(argv)
    shell.main(shell.getArgs())
    return shell.getCalcedParams()


if __name__ == "__main__":
    main()

    # CLI examples:
    # Binomial array, print variables:
    # python array_calculator.py binomial_array -N 6 -f 3e9
    #
    # Dolph-Tschebyscheff array with -26 dB sidelobes:
    # python array_calculator.py dolph_tschebyscheff -N 10 -sll 26 --norm edge
    #
    # Triangular taper, save a polar plot to a file (no display needed):
    # python array_calculator.py triangular_array -N 16 --save tri.png --plot-style polar
    #
    # Hamming taper, rectangular plot in a window:
    # python array_calculator.py hamming_array -N 20 --plot --plot-style rect
    #
    # Export the radiation pattern to CSV:
    # python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --csv pattern.csv
    #
    # Taylor n-bar array with -30 dB near-in sidelobes, nbar = 6:
    # python array_calculator.py taylor_array -N 20 -sll 30 -nbar 6 --verbose
    #
    # Return Variables instead of printing:
    # python array_calculator.py binomial_array -N 6 --variable_return