# ArrayCalculator
in progress AntennaCAT compatible array calculator.

## NOT STABLE, NOT COMPLETE

TODO:
* Add links and references
* Update links and references in the text using the [1] placeholder (not everything is from the Balanis book)
* Update table of contents
* Finish adding the arrays from lecture material
* Upload lecture documents for references
* Add in the images to the README for the output combinations
* Set up tests 



This project is a CLI-based linear antenna array synthesis (sort of, full functionality in progress) tool. 
It is written as an AntennaCAT [ADD LINK] compatible tool, similar to [ADD LINK] [AntennaCalculator]().
This tool designs the excitation coefficients for a multi-element array given a design
goal, then reports the beam characteristics and (optionally) exports the radiation pattern.



The calculator features the following synthesis methods:
* Uniform array (reference case: equal excitation, ~ -13.2 dB sidelobes)
* Binomial array (no sidelobes, amplitudes from Pascal's triangle)
* Dolph-Tschebyscheff array (narrowest beam for a specified, equal sidelobe level)
* Taylor n-bar array (near-in sidelobes held near a design level, decaying tails)
* Triangular taper (~ -26 dB first sidelobe; the square of the uniform pattern)
* Cosine taper (~ -23 dB)
* Cosine-squared taper (~ -31 dB)
* Hann taper (~ -31 dB; mathematically identical to cosine-squared)
* Hamming taper (~ -41 dB; pedestal tuned to cancel the nearest sidelobe)
* Blackman taper (~ -58 dB; three-term cosine, lowest sidelobes of the set)


Supported outputs:
* Excitation coefficients (edge- or center-normalized), HPBW, directivity
* Radiation-pattern export to CSV
* Radiation-pattern plot (optional, requires matplotlib): polar, rectangular,
  or both; shown in a window or saved to a file with `--save`



The synthesis methods are based on the standard formulations in [1]and [2]. The
Dolph-Tschebyscheff coefficients are verified numerically to place the pattern
sidelobes at the requested level.



## Table of Contents
* [Requirements](#requirements)
* [File Structure](#file-structure)
* [Usage](#usage)
    * [Binomial Array Usage](#binomial-array-usage)
    * [Dolph-Tschebyscheff Usage](#dolph-tschebyscheff-usage)
* [Examples](#examples)
* [Tests](#tests)
* [References](#references)

## Requirements

This project requires numpy and pint. matplotlib is optional and only needed
for the `--plot` flag; everything else runs without it.

Tested on Python 3.12.



Use 'pip install -r requirements.txt' to install the following dependencies:
```python
numpy
pint
matplotlib   # optional, for --plot only
```

For development, also install:

```python 
pytest
```



## Building and Installing

The project builds with [uv](https://docs.astral.sh/uv/). To build a local
wheel and install it into your environment:

```python
# install uv
pip install uv
# navigate to the ArrayCalculator directory
cd .\ArrayCalculator
# build the package
# a 'dist' directory should be created in ArrayCalculator
uv build
# install the package locally (matches whatever version was just built)
uv pip install dist/arraycalculator-*.whl
```

This installs the array modules so they can be imported directly (e.g.
`from binomial_array import BinomialArray`). The calculator is still run as a
script -- `python src/array_calculator.py binomial_array -N 6 -f 3e9` -- since
no console entry point is configured.



## File Structure

```python

ArrayCalculator/
├── .python-version
├── pyproject.toml
├── README.md
├── LICENSE
├── .gitignore
├── media/
├── src/
│   ├── __init__.py
│   ├── array_calculator.py
│   ├── array_common.py
│   ├── uniform_array.py
│   ├── binomial_array.py
│   ├── dolph_tschebyscheff.py
│   ├── taylor_array.py
│   ├── triangular_array.py
│   ├── cosine_array.py
│   ├── cosine_squared_array.py
│   ├── hann_array.py
│   ├── hamming_array.py
│   ├── blackman_array.py
│   └── windowed_array.py     (staged: consolidated taper module, not yet wired in)
└── tests/
    ├── __init__.py
    ├── test_array_calculator.py (OLD: superseded by the per-array files below)
    ├── test_helpers.py
    ├── test_uniform_array.py
    ├── test_binomial_array.py
    ├── test_dolph_tschebyscheff.py
    ├── test_triangular_array.py
    ├── test_cosine_array.py
    ├── test_cosine_squared_array.py
    ├── test_hann_array.py
    ├── test_hamming_array.py
    └── test_blackman_array.py
```

```
array_calculator.py       Main entry point. Builds the CLI and dispatches
                          to the synthesis classes (ArrayCalculator class).
array_common.py           Shared base class (ArrayCommon): printing helpers,
                          frequency/spacing reporting, array-factor
                          computation, CSV export, and plotting (polar/rect).
uniform_array.py          UniformArray class. Equal excitation reference case.
binomial_array.py         BinomialArray class. Pascal's-triangle amplitudes,
                          HPBW, and directivity.
dolph_tschebyscheff.py    DolphTschebyscheff class. Coefficient synthesis for
                          a specified sidelobe level, with edge/center
                          normalization.
taylor_array.py           TaylorArray class. n-bar line-source distribution.
triangular_array.py       TriangularArray class. Linear taper, ~ -26 dB.
cosine_array.py           CosineArray class. Half-cosine taper, ~ -23 dB.
cosine_squared_array.py   CosineSquaredArray class. Raised cosine, ~ -31 dB.
hann_array.py             HannArray class. Identical to cosine-squared.
hamming_array.py          HammingArray class. 0.54/0.46 raised cosine, ~ -41 dB.
blackman_array.py         BlackmanArray class. Three-term cosine, ~ -58 dB.
windowed_array.py         Staged WindowedArray class: collapses the six
                          closed-form tapers into one taper table. Not yet
                          imported by the CLI; the standalone files remain
                          the source of truth.
```

All array classes inherit from `ArrayCommon`, so each method only supplies its
excitation amplitudes; the shared pattern, output, and plotting machinery lives
in one place.

## Usage

```
usage: array_calculator.py [--help] [--version]
                           {uniform_array,binomial_array,triangular_array,
                            cosine_array,cosine_squared_array,hann_array,
                            hamming_array,blackman_array,dolph_tschebyscheff,
                            taylor_array}
                           ...

Array Calculator

positional arguments:
  {uniform_array, binomial_array, triangular_array, cosine_array,
   cosine_squared_array, hann_array, hamming_array, blackman_array,
   dolph_tschebyscheff, taylor_array}
                        sub-command help

optional arguments:
  --help                Show this help message and exit
  --version             show program's version number and exit
```

### Common arguments

Every subcommand accepts the same base arguments below. The two parametric
methods (`dolph_tschebyscheff`, `taylor_array`) add a few of their own on top
of these.

```
  --help                Show this help message and exit
  --verbose             Print extra synthesis detail (raw amplitudes, etc.)
  -N, --elements N      Number of array elements (required)
  -f, --frequency F     Frequency in Hz (optional; gives physical spacing)
  -d, --spacing D       Element spacing as a fraction of wavelength (default 0.5)
  --scan SCAN           Beam steering angle in degrees (default 90 = broadside)
  -u, --unit {meter,centimeter,millimeter,inch}
                        Unit of measurement (default centimeter)
  --csv CSV             Export radiation pattern to CSV file
  --plot                Plot radiation pattern (requires matplotlib)
  --plot-style {polar,rect,both}
                        Plot layout: polar dial, rectangular dB, or both
                        (default both)
  --save SAVE           Save the plot to this path instead of opening a window
                        (implies --plot; works without a display)
  --variable_return     Return variables instead of printing
```

### Closed-form taper arrays

`uniform_array`, `binomial_array`, `triangular_array`, `cosine_array`,
`cosine_squared_array`, `hann_array`, `hamming_array`, and `blackman_array` take
only the common arguments. For example:

```
python array_calculator.py triangular_array -N 16
python array_calculator.py blackman_array -N 20 --plot --plot-style polar
python array_calculator.py hamming_array -N 16 --save hamming.png
```

### Dolph-Tschebyscheff Usage

Adds `-sll` (required) and `--norm` to the common arguments.

```
usage: array_calculator.py dolph_tschebyscheff [common args]
                                               -sll SIDELOBE_LEVEL
                                               [--norm {edge,center}]

  -sll, --sidelobe_level SIDELOBE_LEVEL
                        Desired sidelobe level in dB (e.g. 26, 30) (required)
  --norm {edge,center}  Normalize to edge or center element (default edge)
```

### Taylor Array Usage

Adds `-sll` (required) and `-nbar` to the common arguments.

```
usage: array_calculator.py taylor_array [common args]
                                        -sll SIDELOBE_LEVEL [-nbar NBAR]

  -sll, --sidelobe_level SIDELOBE_LEVEL
                        Desired (near-in) sidelobe level in dB (required)
  -nbar, --nbar NBAR    Number of near-in sidelobes held near the design level
                        (default 5)
```

## Examples

Each synthesis method is its own subcommand. The examples below show
representative output; element amplitudes, HPBW, and directivity are printed
for every array, with method-specific extras (sidelobe level, voltage ratio,
etc.) where they apply. Add `--verbose` to any command for the raw,
un-normalized amplitudes and intermediate quantities.

**Uniform array (reference case), 10 elements at 3 GHz:**
```
python array_calculator.py uniform_array -N 10 -f 3e9

[*] N = 10
[*] Amplitudes (normalized) = [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000]
[*] Wavelength = 10.00 centimeter
[*] Element spacing d = 5.00 centimeter
[*] HPBW = 10.21 deg
[*] Directivity = 10
[*] Directivity = 10.00 dB
```

At d = lambda/2 the uniform array's directivity is simply N, and its first
sidelobe sits at the classic -13.2 dB -- the narrowest beam and highest
sidelobes of the distributions here, making it the natural baseline.

**Binomial array, 6 elements at 3 GHz:**
```
python array_calculator.py binomial_array -N 6 -f 3e9 --verbose

[*] N = 6
[*] Amplitudes (raw) = [1.000, 5.000, 10.000, 10.000, 5.000, 1.000]
[*] Amplitudes (normalized) = [0.100, 0.500, 1.000, 1.000, 0.500, 0.100]
[*] Wavelength = 10.00 centimeter
[*] Element spacing d = 5.00 centimeter
[*] HPBW = 27.16 deg
[*] Directivity = 4.34
[*] Directivity = 6.37 dB
```

**Dolph-Tschebyscheff, 10 elements, -26 dB sidelobes:**
```
python array_calculator.py dolph_tschebyscheff -N 10 -sll 26 --norm edge

[*] N = 10
[*] Sidelobe level = 26.00 dB
[*] Amplitudes (edge-normalized) = [1.000, 1.355, 1.968, 2.479, 2.769, 2.769, 2.479, 1.968, 1.355, 1.000]
[*] Element spacing d = 0.50 lambda
[*] Directivity = 18.16
[*] Directivity = 12.59 dB
```



**Triangular taper, 10 elements at 3 GHz:**
```
python array_calculator.py triangular_array -N 10 -f 3e9

[*] N = 10
[*] Amplitudes (normalized) = [0.200, 0.400, 0.600, 0.800, 1.000, 1.000, 0.800, 0.600, 0.400, 0.200]
[*] Wavelength = 10.00 centimeter
[*] Element spacing d = 5.00 centimeter
[*] HPBW = 13.46 deg
[*] Directivity = 8.18
[*] Directivity = 9.13 dB
```

**Taylor n-bar array, 20 elements, -30 dB near-in sidelobes, nbar = 6:**
```
python array_calculator.py taylor_array -N 20 -sll 30 -nbar 6 --verbose

[*] N = 20
[*] Sidelobe level = 30.00 dB
[*] nbar = 6
[*] R (voltage ratio) = 31.62
[*] sigma = 1.06
[*] Amplitudes (raw) = [0.414, 0.468, 0.586, 0.757, 0.945, 1.119, 1.273, 1.403, 1.495, 1.541, 1.541, 1.495, 1.403, 1.273, 1.119, 0.945, 0.757, 0.586, 0.468, 0.414]
[*] Amplitudes (normalized) = [0.269, 0.304, 0.380, 0.492, 0.613, 0.726, 0.826, 0.910, 0.970, 1.000, 1.000, 0.970, 0.910, 0.826, 0.726, 0.613, 0.492, 0.380, 0.304, 0.269]
[*] Element spacing d = 0.50 lambda
[*] HPBW = 6.42 deg
[*] Directivity = 17.17
[*] Directivity = 12.35 dB
```

Unlike Dolph-Tschebyscheff, whose sidelobes are all equal, the Taylor
distribution holds only the first `nbar` sidelobes near the design level and
lets the rest decay -- trading a slightly wider beam for better aperture
efficiency.

**Hamming taper, 16 elements:**
```
python array_calculator.py hamming_array -N 16

[*] N = 16
[*] Amplitudes (normalized) = [0.090, 0.159, 0.287, 0.454, 0.635, 0.803, 0.931, 1.000, 1.000, 0.931, 0.803, 0.635, 0.454, 0.287, 0.159, 0.090]
[*] Element spacing d = 0.50 lambda
[*] HPBW = 9.35 deg
[*] Directivity = 11.74
[*] Directivity = 10.70 dB
```

**Blackman taper, 16 elements:**
```
python array_calculator.py blackman_array -N 16

[*] N = 16
[*] Amplitudes (normalized) = [0.004, 0.035, 0.113, 0.253, 0.451, 0.678, 0.880, 1.000, 1.000, 0.880, 0.678, 0.451, 0.253, 0.113, 0.035, 0.004]
[*] Element spacing d = 0.50 lambda
[*] HPBW = 11.80 deg
[*] Directivity = 9.27
[*] Directivity = 9.67 dB
```

The remaining tapers -- `cosine_array`, `cosine_squared_array`, and
`hann_array` -- take the same arguments and print output in the same form.
(Cosine-squared and Hann are the same distribution and produce identical
results; see their module docstrings.)

**Return variables instead of printing** (for use as a library):
```
python array_calculator.py binomial_array -N 6 --variable_return
```

With `--variable_return` the calculator suppresses printing and the
`ArrayCalculator` object exposes the computed values via `getCalcedParams()`
-- handy when driving the synthesis from other code rather than the terminal.


**Export the pattern to CSV and steer the beam to 60 degrees:**
```
python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --csv pattern.csv --scan 60
```

The CSV contains three columns: `theta_deg`, `AF_linear`, and `AF_dB`:

```
theta_deg,AF_linear,AF_dB
0.000,0.000000,-240.0000
0.250,0.000000,-130.8449
0.500,0.000001,-118.8037
...
90.000,1.000000,0.0000
...
```


**Matplotlib visualization:**

The `--plot` flag renders the radiation pattern; `--plot-style` selects the
layout and `--save` writes it to a file (which also works on a headless machine,
where opening a window would otherwise fail):

```
# open an interactive window with both polar and rectangular views
python array_calculator.py dolph_tschebyscheff -N 10 -sll 26 --plot

# just the polar dial
python array_calculator.py hamming_array -N 20 --plot --plot-style polar

# save a rectangular dB plot to a file, no display needed
python array_calculator.py blackman_array -N 24 --save blackman.png --plot-style rect
```

The polar view shows the 0-180 deg half-plane mirrored about boresight (the
beam points up at 90 deg); the rectangular view plots normalized |AF| in dB
against theta. When a sidelobe level applies (Dolph, Taylor), it is drawn as a
dashed reference line.

<!-- Add output images below. Suggested captions: -->

<!-- Dolph-Tschebyscheff, both views -->
![Dolph-Tschebyscheff radiation pattern, polar and rectangular](media/dolph_both.png)

<!-- A low-sidelobe taper, polar -->
![Blackman radiation pattern](media/blackman_polar.png)


## Tests

The test suite is in `test_array_calculator.py` and uses `pytest`. It checks the
synthesis math against the worked examples in [1] and verifies the computed
patterns physically (e.g. that a Dolph-Tschebyscheff design's sidelobes actually
land at the requested level, and that a binomial array has none).

The suite is split one file per array (`test_uniform_array.py`,
`test_binomial_array.py`, `test_triangular_array.py`, and so on), with shared
helpers in `test_helpers.py`. The older monolithic `test_array_calculator.py` is
kept for reference but superseded by these.

Run the full suite with:

```
pytest -v
```

Or a single array's tests directly:

```
python test_triangular_array.py
```

Coverage includes:
* Binomial amplitudes (Pascal's triangle), HPBW, and directivity
* Dolph-Tschebyscheff z0, R, and coefficient synthesis for even and odd N
* Edge- vs. center-normalization consistency
* Closed-form tapers: symmetry, center peak, and the characteristic peak
  sidelobe level of each (triangular ~ -26 dB, cosine ~ -23 dB,
  cosine-squared / Hann ~ -31 dB, Hamming ~ -41 dB, Blackman ~ -58 dB)
* The Hann / cosine-squared equivalence (identical amplitudes across N)
* Physical pattern checks: target sidelobe level, equiripple sidelobes,
  no-sidelobe binomial behavior, and beam steering
* CSV export format
* CLI dispatch, defaults, `--variable_return`, and argument validation

Running `pytest` requires `pytest` in addition to the project's `numpy` and
`pint` dependencies.


## References

[1]: C. A. Balanis, Antenna Theory: Analysis and Design. Hoboken, New Jersey: Wiley, 2016.
[2]: 
[3]: