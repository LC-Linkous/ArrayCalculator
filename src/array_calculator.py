#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  array_calculator.py
#
#   This is the main entry point of the array calculator.
#   This project is designed as an AntennaCAT compatible calcuator.
#   Refer to the README for argument formats.
##--------------------------------------------------------------------\


import argparse

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


def _add_common_args(p, require_elements=True):
    # The argument set every array subcommand shares. Centralized here so
    # all arrays stay in parity (spacing, steering, units, output, plotting)
    # and new flags only have to be added in one place.
    p.add_argument('--help', action='help', default=argparse.SUPPRESS,
                   help='Show this help message and exit')
    p.add_argument('--verbose', action='store_true')
    p.add_argument('-N', '--elements', type=int, required=require_elements,
                   help='Number of array elements')
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
    p.add_argument('--csv', type=str, required=False,
                   help='Export radiation pattern to CSV file')
    p.add_argument('--plot', action='store_true', required=False, default=False,
                   help='Plot radiation pattern (requires matplotlib)')
    p.add_argument('--plot-style', dest='plot_style', type=str,
                   choices=['polar', 'rect', 'both'], required=False, default='both',
                   help='Plot layout: polar dial, rectangular dB, or both (default both)')
    p.add_argument('--save', type=str, required=False, default=None,
                   help='Save the plot to this path instead of opening a window '
                        '(implies --plot; works without a display)')
    p.add_argument('--variable_return', action='store_true', required=False, default=False,
                   help='Return Variables instead of printing')
    p.add_argument('--norm', type=str, choices=['edge', 'center'],
                   required=False, default='center',
                   help='Normalize amplitudes to the edge or center element '
                        '(default center; Dolph-Tschebyscheff conventionally uses edge)')


class ArrayCalculator:
    def __init__(self, a=None):
        main_parser = argparse.ArgumentParser(description='Array Calculator', add_help=False)
        main_parser.add_argument('--help', action='help', default=argparse.SUPPRESS,
                                 help='Show this help message and exit')
        main_parser.add_argument('--version', action='version', version='%(prog)s 1.0')

        subparsers = main_parser.add_subparsers(help='sub-command help',
                                                dest='subparser_name')

        # Arrays that take only the common arguments.
        for name in SIMPLE_ARRAYS:
            sp = subparsers.add_parser(name, add_help=False)
            _add_common_args(sp)

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
        _add_common_args(ev, require_elements=False)
        ev.add_argument('-g', '--geometry', type=str, required=True,
                        help='Geometry CSV: position_lambda[, amplitude][, phase_deg]')

        self.args = main_parser.parse_args(a)
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


if __name__ == "__main__":

    shell = ArrayCalculator()
    args = shell.getArgs()
    shell.main(args)

    
    # ==================================================================
    # CLI EXAMPLES  (each shown two ways: as a shell command, and as the
    # equivalent scripted ArrayCalculator call for automation/screenshots)
    # ------------------------------------------------------------------
    # The scripted form mirrors the __main__ block: build ArrayCalculator
    # with an argv list, then call .main(.getArgs()). Add '--variable_return'
    # to suppress printing and collect the results via .getCalcedParams().
    #
    # NOTE on return shapes (what getCalcedParams() gives back):
    #   most arrays         -> (amps, hpbw_deg, directivity_db)
    #   dolph_tschebyscheff -> (amps, R, z0, directivity_db)      # 4-tuple
    #   woodward_lawson     -> (amps, phase_deg, directivity_db)  # phase, not HPBW
    #   evaluate            -> a results dict (not a tuple)
    # ==================================================================
    

    # --- Uniform array (reference case), physical spacing from frequency ---
    # python array_calculator.py uniform_array -N 10 -f 3e9
    # shell = ArrayCalculator(['uniform_array', '-N', '10', '-f', '3e9'])
    # shell.main(shell.getArgs())
    

    # --- Binomial array, verbose (also prints raw Pascal-triangle amplitudes) ---
    # python array_calculator.py binomial_array -N 6 -f 3e9 --verbose
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '-f', '3e9', '--verbose'])
    # shell.main(shell.getArgs())
    

    # --- Dolph-Tschebyscheff, -26 dB sidelobes, edge-normalized ---
    # python array_calculator.py dolph_tschebyscheff -N 10 -sll 26 --norm edge
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '10', '-sll', '26', '--norm', 'edge'])
    # shell.main(shell.getArgs())
    

    # --- Dolph-Tschebyscheff, center-normalized variant ---
    # python array_calculator.py dolph_tschebyscheff -N 12 -sll 30 --norm center
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '12', '-sll', '30', '--norm', 'center'])
    # shell.main(shell.getArgs())
    

    # --- Taylor n-bar, -30 dB near-in sidelobes, nbar = 6, verbose ---
    # python array_calculator.py taylor_array -N 20 -sll 30 -nbar 6 --verbose
    # shell = ArrayCalculator(['taylor_array', '-N', '20', '-sll', '30', '-nbar', '6', '--verbose'])
    # shell.main(shell.getArgs())
    

    # --- Villeneuve n-bar (discrete Taylor), same design point ---
    # python array_calculator.py villeneuve_array -N 20 -sll 30 -nbar 6 --verbose
    # shell = ArrayCalculator(['villeneuve_array', '-N', '20', '-sll', '30', '-nbar', '6', '--verbose'])
    # shell.main(shell.getArgs())
    

    # --- Triangular taper, save a polar plot to a file (headless-safe) ---
    # python array_calculator.py triangular_array -N 16 --save tri.png --plot-style polar
    # shell = ArrayCalculator(['triangular_array', '-N', '16', '--save', 'tri.png', '--plot-style', 'polar'])
    # shell.main(shell.getArgs())


    # --- Bartlett taper (zero-endpoint triangular) ---
    # python array_calculator.py bartlett_array -N 16
    # shell = ArrayCalculator(['bartlett_array', '-N', '16'])
    # shell.main(shell.getArgs())
    

    # --- Cosine taper ---
    # python array_calculator.py cosine_array -N 12
    # shell = ArrayCalculator(['cosine_array', '-N', '12'])
    # shell.main(shell.getArgs())
    

    # --- Hann taper ---
    # python array_calculator.py hann_array -N 16
    # shell = ArrayCalculator(['hann_array', '-N', '16'])
    # shell.main(shell.getArgs())
    

    # --- Hamming taper, rectangular plot in an interactive window ---
    # python array_calculator.py hamming_array -N 20 --plot --plot-style rect
    # shell = ArrayCalculator(['hamming_array', '-N', '20', '--plot', '--plot-style', 'rect'])
    # shell.main(shell.getArgs())
    

    # --- Blackman taper, both plot views saved to a file ---
    # python array_calculator.py blackman_array -N 24 --save blackman.png
    # shell = ArrayCalculator(['blackman_array', '-N', '24', '--save', 'blackman.png'])
    # shell.main(shell.getArgs())
    

    # --- Kaiser taper, beta = 8 (low sidelobes) ---
    # python array_calculator.py kaiser_array -N 16 -beta 8
    # shell = ArrayCalculator(['kaiser_array', '-N', '16', '-beta', '8'])
    # shell.main(shell.getArgs())
    

    # --- Kaiser taper, beta = 4 (closer to uniform) ---
    # python array_calculator.py kaiser_array -N 20 -beta 4
    # shell = ArrayCalculator(['kaiser_array', '-N', '20', '-beta', '4'])
    # shell.main(shell.getArgs())
    
    
    # --- Woodward-Lawson flat-top beam over a +/-30 deg sector ---
    # python array_calculator.py woodward_lawson_array -N 20 --shape flat_top --sector 30
    # shell = ArrayCalculator(['woodward_lawson_array', '-N', '20', '--shape', 'flat_top', '--sector', '30'])
    # shell.main(shell.getArgs())
    

    # --- Woodward-Lawson cosecant-squared coverage beam, verbose (prints phases) ---
    # python array_calculator.py woodward_lawson_array -N 24 --shape cosecant_squared --sector 35 --verbose
    # shell = ArrayCalculator(['woodward_lawson_array', '-N', '24', '--shape', 'cosecant_squared', '--sector', '35', '--verbose'])
    # shell.main(shell.getArgs())
    

    # --- Steer the beam to 60 degrees (works on any array) ---
    # python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --scan 60 --verbose
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '8', '-sll', '25', '--scan', '60', '--verbose'])
    # shell.main(shell.getArgs())
    

    # --- Ordinary endfire: steer to the array axis (use closer spacing) ---
    # python array_calculator.py uniform_array -N 10 -d 0.25 --scan 0
    # shell = ArrayCalculator(['uniform_array', '-N', '10', '-d', '0.25', '--scan', '0'])
    # shell.main(shell.getArgs())
    

    # --- Export the radiation pattern to CSV (with steering) ---
    # python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --csv pattern.csv --scan 60
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '8', '-sll', '25', '--csv', 'pattern.csv', '--scan', '60'])
    # shell.main(shell.getArgs())
    

    # --- Evaluate an arbitrary geometry from a CSV (the optimizer hook) ---
    # python array_calculator.py evaluate -g ./example_data/geometry.csv
    # shell = ArrayCalculator(['evaluate', '-g', './example_data/geometry.csv'])
    # shell.main(shell.getArgs())
    

    # --- Return variables instead of printing, then read them back ---
    # python array_calculator.py binomial_array -N 6 --variable_return
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '--variable_return'])
    # shell.main(shell.getArgs())
    # amps, hpbw_deg, directivity_db = shell.getCalcedParams()
    # print("done!")
    

    #   # Dolph returns a 4-tuple (note the extra R, z0):
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '10', '-sll', '26', '--variable_return'])
    # shell.main(shell.getArgs())
    # amps, R, z0, directivity_db = shell.getCalcedParams()
    # print("done!")
    
    #   # evaluate returns a results dict, not a tuple:
    # shell = ArrayCalculator(['evaluate', '-g', './example_data/geometry.csv', '--variable_return'])
    # shell.main(shell.getArgs())
    # results = shell.getCalcedParams()        # results['peak_sidelobe_db'], etc.
    # print(results)
