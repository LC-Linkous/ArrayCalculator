#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  array_calculator.py
#
#   This is the main entry point of the array calculator. It mirrors the
#   structure of antenna_calculator.py so the two can be used and imported
#   side by side. Refer to the README for argument formats.
##--------------------------------------------------------------------\


import argparse
from binomial_array import BinomialArray
from dolph_tschebyscheff import DolphTschebyscheff


class ArrayCalculator:
    def __init__(self, a=None):
        main_parser = argparse.ArgumentParser(description='Array Calculator', add_help=False)
        main_parser.add_argument('--help', action='help', default=argparse.SUPPRESS,
                                 help='Show this help message and exit')
        main_parser.add_argument('--version', action='version', version='%(prog)s 1.0')

        subparsers = main_parser.add_subparsers(help='sub-command help', dest='subparser_'
                                                                              'name')

        # BINOMIAL ARRAY
        binomial_array_subparser = subparsers.add_parser('binomial_array', add_help=False)
        binomial_array_subparser.add_argument('--help', action='help', default=argparse.SUPPRESS,
                                              help='Show this help message and exit')
        binomial_array_subparser.add_argument('--verbose', action='store_true')
        binomial_array_subparser.add_argument('-N', '--elements', type=int, required=True,
                                              help='Number of array elements')
        binomial_array_subparser.add_argument('-f', '--frequency', type=float, required=False,
                                              help='Frequency in Hz (optional; gives physical spacing)')
        binomial_array_subparser.add_argument('-d', '--spacing', type=float, required=False, default=0.5,
                                              help='Element spacing as a fraction of wavelength (default 0.5)')
        binomial_array_subparser.add_argument('--scan', type=float, required=False, default=90.0,
                                              help='Beam steering angle in degrees (default 90 = broadside)')
        binomial_array_subparser.add_argument('-u', '--unit', type=str,
                                              choices=['meter', 'centimeter', 'millimeter', 'inch'],
                                              required=False, default='centimeter',
                                              help='Unit of measurement')
        binomial_array_subparser.add_argument('--csv', type=str, required=False,
                                              help='Export radiation pattern to CSV file')
        binomial_array_subparser.add_argument('--plot', action='store_true', required=False, default=False,
                                              help='Plot radiation pattern (requires matplotlib)')
        binomial_array_subparser.add_argument('--variable_return', action='store_true', required=False, default=False,
                                              help='Return Variables instead of printing')

        # DOLPH-TSCHEBYSCHEFF ARRAY
        dolph_tschebyscheff_subparser = subparsers.add_parser('dolph_tschebyscheff', add_help=False)
        dolph_tschebyscheff_subparser.add_argument('--help', action='help', default=argparse.SUPPRESS,
                                                   help='Show this help message and exit')
        dolph_tschebyscheff_subparser.add_argument('--verbose', action='store_true')
        dolph_tschebyscheff_subparser.add_argument('-N', '--elements', type=int, required=True,
                                                   help='Number of array elements')
        dolph_tschebyscheff_subparser.add_argument('-sll', '--sidelobe_level', type=float, required=True,
                                                   help='Desired sidelobe level in dB (e.g. 26, 30)')
        dolph_tschebyscheff_subparser.add_argument('--norm', type=str, choices=['edge', 'center'],
                                                   required=False, default='edge',
                                                   help='Normalize to edge or center element (default edge)')
        dolph_tschebyscheff_subparser.add_argument('-f', '--frequency', type=float, required=False,
                                                   help='Frequency in Hz (optional; gives physical spacing)')
        dolph_tschebyscheff_subparser.add_argument('-d', '--spacing', type=float, required=False, default=0.5,
                                                   help='Element spacing as a fraction of wavelength (default 0.5)')
        dolph_tschebyscheff_subparser.add_argument('--scan', type=float, required=False, default=90.0,
                                                   help='Beam steering angle in degrees (default 90 = broadside)')
        dolph_tschebyscheff_subparser.add_argument('-u', '--unit', type=str,
                                                   choices=['meter', 'centimeter', 'millimeter', 'inch'],
                                                   required=False, default='centimeter',
                                                   help='Unit of measurement')
        dolph_tschebyscheff_subparser.add_argument('--csv', type=str, required=False,
                                                   help='Export radiation pattern to CSV file')
        dolph_tschebyscheff_subparser.add_argument('--plot', action='store_true', required=False, default=False,
                                                   help='Plot radiation pattern (requires matplotlib)')
        dolph_tschebyscheff_subparser.add_argument('--variable_return', action='store_true', required=False, default=False,
                                                   help='Return Variables instead of printing')

        self.args = main_parser.parse_args(a)
        self.calcedParams = None  # to catch returned vars if they exist

    def main(self, args):
        if args.subparser_name == 'binomial_array':
            b = BinomialArray(args)
            self.calcedParams = b.binomial_array_calculator()

        if args.subparser_name == 'dolph_tschebyscheff':
            dt = DolphTschebyscheff(args)
            self.calcedParams = dt.dolph_tschebyscheff_calculator()

    def getArgs(self):
        return self.args

    def getCalcedParams(self):
        return self.calcedParams


if __name__ == "__main__":

    shell = ArrayCalculator()
    args = shell.getArgs()
    shell.main(args)

    # CLI examples:
    # Binomial array, print variables:
    # python array_calculator.py binomial_array -N 6 -f 3e9
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '-f', '3e9'])
    # shell.main(shell.getArgs())

    # Binomial array, print all variables:
    # python array_calculator.py binomial_array -N 6 -f 3e9 --verbose
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '-f', '3e9', '--verbose'])
    # shell.main(shell.getArgs())

    # Dolph-Tschebyscheff array with -26 dB sidelobes:
    # python array_calculator.py dolph_tschebyscheff -N 10 -sll 26 --norm edge
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '10', '-sll','26', '--norm', 'edge'])
    # shell.main(shell.getArgs())

    # Export the radiation pattern to CSV:
    # python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --csv pattern.csv
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '8', '-sll','25', '--csv', 'pattern.csv'])
    # shell.main(shell.getArgs())

    # Return Variables instead of printing:
    # python array_calculator.py binomial_array -N 6 --variable_return
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '--variable_return'])
    # shell.main(shell.getArgs())