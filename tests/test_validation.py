#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_validation.py
#
#   Bad CLI inputs must stop with a clear argparse error (SystemExit 2)
#   rather than crash with ZeroDivisionError / IndexError, print NaN, or
#   silently produce a corrupted array. The evaluator must reject the
#   synthesis-only flags it would otherwise ignore.
##--------------------------------------------------------------------\

import numpy as np
import pytest

from array_calculator import ArrayCalculator
from villeneuve_array import VilleneuveArray
from dolph_tschebyscheff import DolphTschebyscheff
from test_helpers import make_args


BAD = [
    (["binomial_array", "-N", "1"], "at least 2"),
    (["dolph_tschebyscheff", "-N", "1", "-sll", "26"], "at least 2"),
    (["uniform_array", "-N", "8", "-d", "0"], "spacing"),
    (["uniform_array", "-N", "8", "-d", "-0.5"], "spacing"),
    (["uniform_array", "-N", "8", "-f", "0"], "frequency"),
    (["uniform_array", "-N", "8", "--scan", "200"], "scan"),
    (["dolph_tschebyscheff", "-N", "8", "-sll", "0"], "-sll"),
    (["dolph_tschebyscheff", "-N", "8", "-sll", "-20"], "-sll"),
    (["taylor_array", "-N", "16", "-sll", "30", "-nbar", "0"], "nbar"),
    (["villeneuve_array", "-N", "10", "-sll", "30", "-nbar", "6"], "nbar"),
    (["villeneuve_array", "-N", "10", "-sll", "30", "-nbar", "12"], "nbar"),
    (["kaiser_array", "-N", "16", "-beta", "-1"], "beta"),
    (["woodward_lawson_array", "-N", "16", "--sector", "0"], "sector"),
    (["woodward_lawson_array", "-N", "16", "--floor", "-0.1"], "floor"),
    (["woodward_lawson_array", "-N", "16", "--shape", "cosecant_squared"], "scan"),
]


@pytest.mark.parametrize("argv,needle", BAD)
def test_bad_input_is_a_clean_cli_error(argv, needle, capsys):
    with pytest.raises(SystemExit) as e:
        ArrayCalculator(argv)
    assert e.value.code == 2
    assert needle in capsys.readouterr().err


@pytest.mark.parametrize("flag", [["-N", "8"], ["-d", "0.7"], ["--scan", "60"],
                                  ["-f", "3e9"], ["--norm", "edge"]])
def test_evaluate_rejects_ignored_flags(flag, tmp_path):
    g = tmp_path / "g.csv"
    g.write_text("position_lambda\n0\n0.5\n1.0\n")
    with pytest.raises(SystemExit):
        ArrayCalculator(["evaluate", "-g", str(g)] + flag)


def test_no_subcommand_is_an_error():
    with pytest.raises(SystemExit):
        ArrayCalculator([])


def test_good_edge_values_accepted():
    # Boundaries that are legal must still run.
    for argv in (["uniform_array", "-N", "2"],
                 ["villeneuve_array", "-N", "10", "-sll", "30", "-nbar", "5"],
                 ["villeneuve_array", "-N", "11", "-sll", "30", "-nbar", "6"],
                 ["kaiser_array", "-N", "8", "-beta", "0"],
                 ["uniform_array", "-N", "8", "--scan", "0"],
                 ["woodward_lawson_array", "-N", "16", "--shape",
                  "cosecant_squared", "--scan", "80"]):
        shell = ArrayCalculator(argv + ["--variable_return"])
        shell.main(shell.getArgs())
        assert shell.getCalcedParams() is not None


def test_programmatic_guards():
    v = VilleneuveArray(make_args(elements=10))
    with pytest.raises(ValueError):
        v.amplitudes(10, 30, 6)
    d = DolphTschebyscheff(make_args(elements=1))
    with pytest.raises(ValueError):
        d.amplitudes(1, 26)
    with pytest.raises(ValueError):
        d.amplitudes(8, 0)


def test_villeneuve_max_nbar_still_hits_sll():
    v = VilleneuveArray(make_args(elements=11))
    amps, _, _ = v.amplitudes(11, 25, 6)
    theta, _, af_db = v.pattern_sweep(amps, n_points=20001)
    assert v.peak_sidelobe(theta, af_db) < -24.0
