#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  test_plots.py
#
#   Exercises the plotting path for every array across every plot style.
#   For each (array, style) combination it drives the real CLI with
#   --save (the headless-safe path; --plot's interactive plt.show() can't
#   run without a display), then checks that a valid, non-empty PNG was
#   produced. Run it to regenerate and eyeball the full matrix of plots:
#
#       pytest -v test_plots.py
#       python test_plots.py            # same, plus a printed summary
#
#   Saved images land in a temp dir under pytest; running as __main__
#   writes them to ./plot_out/ so you can open them.
##--------------------------------------------------------------------\

import os

import matplotlib
matplotlib.use("Agg")  # headless backend: no display needed, just files

import pytest

from array_calculator import ArrayCalculator


# Each array and the extra CLI args it needs beyond -N. The parametric
# methods require -sll (and Taylor takes -nbar); the closed-form tapers
# take nothing extra.
ARRAY_ARGS = {
    "uniform_array": [],
    "binomial_array": [],
    "triangular_array": [],
    "cosine_array": [],
    "cosine_squared_array": [],
    "hann_array": [],
    "hamming_array": [],
    "blackman_array": [],
    "dolph_tschebyscheff": ["-sll", "26"],
    "taylor_array": ["-sll", "30", "-nbar", "5"],
}

PLOT_STYLES = ["polar", "rect", "both"]

N = "16"


def _run_one(array_name, style, out_dir):
    """Drive the CLI to save one (array, style) plot; return the path."""
    out_path = os.path.join(out_dir, "%s_%s.png" % (array_name, style))
    argv = [array_name, "-N", N] + ARRAY_ARGS[array_name] + \
           ["--save", out_path, "--plot-style", style]
    shell = ArrayCalculator(argv)
    shell.main(shell.getArgs())
    return out_path


def _is_valid_png(path):
    """Exists, non-trivial size, and starts with the PNG magic number."""
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        return False
    with open(path, "rb") as f:
        return f.read(8) == b"\x89PNG\r\n\x1a\n"


# Full matrix: every array x every plot style.
@pytest.mark.parametrize("array_name", list(ARRAY_ARGS))
@pytest.mark.parametrize("style", PLOT_STYLES)
def test_plot_combination_renders(array_name, style, tmp_path):
    out = _run_one(array_name, style, str(tmp_path))
    assert _is_valid_png(out), \
        "no valid PNG for %s / %s" % (array_name, style)


def test_save_implies_plot(tmp_path):
    # --save alone (no --plot) should still produce a file.
    out = os.path.join(str(tmp_path), "implies.png")
    shell = ArrayCalculator(["uniform_array", "-N", N, "--save", out])
    assert shell.getArgs().plot is True
    shell.main(shell.getArgs())
    assert _is_valid_png(out)


def test_default_plot_style_is_both(tmp_path):
    # No --plot-style given: should default to the combined view.
    shell = ArrayCalculator(["cosine_array", "-N", N,
                             "--save", os.path.join(str(tmp_path), "d.png")])
    assert shell.getArgs().plot_style == "both"


if __name__ == "__main__":
    # Generate the whole matrix into ./plot_out/ and print a summary so you
    # can see every combination at a glance, then open the images.
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_out")
    os.makedirs(out_dir, exist_ok=True)

    print("Generating %d plots (%d arrays x %d styles) into %s\n"
          % (len(ARRAY_ARGS) * len(PLOT_STYLES),
             len(ARRAY_ARGS), len(PLOT_STYLES), out_dir))
    print("%-22s %-6s %-8s %s" % ("array", "style", "size", "file"))
    print("-" * 62)

    ok = 0
    for array_name in ARRAY_ARGS:
        for style in PLOT_STYLES:
            path = _run_one(array_name, style, out_dir)
            valid = _is_valid_png(path)
            ok += valid
            size = "%dk" % (os.path.getsize(path) // 1024) if os.path.exists(path) else "--"
            mark = "OK" if valid else "FAIL"
            print("%-22s %-6s %-8s %s  [%s]"
                  % (array_name, style, size, os.path.basename(path), mark))

    total = len(ARRAY_ARGS) * len(PLOT_STYLES)
    print("\n%d/%d combinations rendered." % (ok, total))