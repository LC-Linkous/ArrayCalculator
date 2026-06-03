#! /usr/bin/python3

##--------------------------------------------------------------------\
#   ArrayCalculator  tests/conftest.py
#
#   pytest auto-imports this file before collecting any tests, so it runs
#   early enough to fix the import path. The array modules live in src/,
#   but the test files import them flat (e.g. `from binomial_array import
#   BinomialArray`). Adding src/ to sys.path here lets those imports
#   resolve no matter what directory pytest is launched from.
#
#   This file belongs in the tests/ directory, as a sibling of the
#   test_*.py files. It expects the project layout:
#
#       ArrayCalculator/
#       ├── src/      <- array modules
#       └── tests/    <- this file + the test_*.py files
#
#   If your array modules live somewhere else relative to tests/, change
#   the SRC path below to match.
##--------------------------------------------------------------------\

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(SRC))