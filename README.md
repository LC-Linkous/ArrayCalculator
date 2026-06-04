# ArrayCalculator
in progress AntennaCAT compatible array calculator.

## NOT STABLE, NOT COMPLETE


This project is a CLI-based linear antenna array synthesis (sort of, full functionality in progress) tool. 
It is written as an [AntennaCAT] (https://github.com/LC-Linkous/AntennaCalculationAutotuningTool) compatible tool, similar to the CLI interface of the [AntennaCalculator](https://github.com/Dollarhyde/AntennaCalculator).
This tool designs the excitation coefficients for a multi-element array given a design
goal, then reports the beam characteristics and (optionally) exports the radiation pattern.



The calculator features the following synthesis methods:
* Uniform array (reference case: equal excitation, ~ -13.2 dB sidelobes)
* Binomial array (no sidelobes, amplitudes from Pascal's triangle)
* Dolph-Tschebyscheff array (narrowest beam for a specified, equal sidelobe level)
* Taylor n-bar array (near-in sidelobes held near a design level, decaying tails)
* Villeneuve n-bar array (discrete-element analog of Taylor; more accurate near-in
  sidelobes for small/moderate N)
* Triangular taper (~ -26 dB first sidelobe; the square of the uniform pattern)
* Bartlett taper (zero-endpoint triangular, ~ -26 dB)
* Cosine taper (~ -23 dB)
* Cosine-squared taper (~ -31 dB)
* Hann taper (~ -31 dB; mathematically identical to cosine-squared)
* Hamming taper (~ -41 dB; pedestal tuned to cancel the nearest sidelobe)
* Blackman taper (~ -58 dB; three-term cosine, lowest sidelobes of the set)
* Kaiser taper (parametric; one shape parameter beta trades beamwidth against
  sidelobe level continuously)
* Woodward-Lawson shaped-beam synthesis (approximates a desired pattern shape --
  flat-top or cosecant-squared -- rather than a pencil beam)

It also provides an **evaluator** that scores an arbitrary array geometry
(element positions, and optionally per-element amplitudes and phases) rather than
synthesizing one, which is the intended hook for an external optimizer.


Supported outputs:
* Excitation coefficients (edge- or center-normalized), HPBW, directivity
* Radiation-pattern export to CSV
* Radiation-pattern plot (optional, requires matplotlib): polar, rectangular,
  or both; shown in a window or saved to a file with `--save`



The synthesis methods are based on the standard formulations (see [references](#references)). The
parametric designs (Dolph-Tschebyscheff, Taylor, Villeneuve) are verified
numerically to place their pattern sidelobes at the requested level, and the
Woodward-Lawson designs are checked to approximate the requested pattern shape.



## Table of Contents
* [Requirements](#requirements)
* [Building and Installing](#building-and-installing)
* [File Structure](#file-structure)
* [Usage](#usage)
    * [Common arguments](#common-arguments)
    * [Beam steering and endfire](#beam-steering-and-endfire)
    * [Closed-form taper arrays](#closed-form-taper-arrays)
    * [Dolph-Tschebyscheff Usage](#dolph-tschebyscheff-usage)
    * [Taylor Array Usage](#taylor-array-usage)
    * [Villeneuve Array Usage](#villeneuve-array-usage)
    * [Kaiser Array Usage](#kaiser-array-usage)
    * [Woodward-Lawson Usage](#woodward-lawson-usage)
    * [Evaluating an arbitrary geometry](#evaluating-an-arbitrary-geometry)
* [Examples](#examples)
* [Tests](#tests)
* [Development and future work](#development-and-future-work)
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
│   ├── array_evaluator.py
│   ├── uniform_array.py
│   ├── binomial_array.py
│   ├── dolph_tschebyscheff.py
│   ├── taylor_array.py
│   ├── villeneuve_array.py
│   ├── triangular_array.py
│   ├── bartlett_array.py
│   ├── cosine_array.py
│   ├── cosine_squared_array.py
│   ├── hann_array.py
│   ├── hamming_array.py
│   ├── blackman_array.py
│   ├── kaiser_array.py
│   ├── woodward_lawson_array.py
│   └── windowed_array.py     (staged: consolidated taper module, not yet wired in)
└── tests/
    ├── __init__.py
    ├── conftest.py            (puts src/ on sys.path for the flat imports)
    ├── test_helpers.py        (shared fixtures only; collects no tests)
    ├── test_array_calculator.py (cross-cutting: array factor, CSV, CLI, edge cases)
    ├── test_array_evaluator.py
    ├── test_normalization.py
    ├── test_uniform_array.py
    ├── test_binomial_array.py
    ├── test_dolph_tschebyscheff.py
    ├── test_taylor_array.py
    ├── test_villeneuve_array.py
    ├── test_triangular_array.py
    ├── test_bartlett_array.py
    ├── test_cosine_array.py
    ├── test_cosine_squared_array.py
    ├── test_hann_array.py
    ├── test_hamming_array.py
    ├── test_kaiser_array.py
    ├── test_woodward_lawson_array.py
    └── test_blackman_array.py
```

```
array_calculator.py       Main entry point. Builds the CLI and dispatches
                          to the synthesis classes (ArrayCalculator class).
array_common.py           Shared base class (ArrayCommon): printing helpers,
                          frequency/spacing reporting, array-factor
                          computation, CSV export, and plotting (polar/rect).
array_evaluator.py        ArrayEvaluator: scores an arbitrary geometry
                          (positions, amplitudes, phases) rather than
                          synthesizing one. The optimizer hook -- exposes
                          evaluate(...) and an `evaluate` CLI subcommand.
uniform_array.py          UniformArray class. Equal excitation reference case.
binomial_array.py         BinomialArray class. Pascal's-triangle amplitudes,
                          HPBW, and directivity.
dolph_tschebyscheff.py    DolphTschebyscheff class. Coefficient synthesis for
                          a specified sidelobe level, with edge/center
                          normalization.
taylor_array.py           TaylorArray class. n-bar line-source distribution.
villeneuve_array.py       VilleneuveArray class. Discrete-element analog of
                          Taylor (stretched Dolph zeros, decaying tails).
triangular_array.py       TriangularArray class. Linear taper, ~ -26 dB.
bartlett_array.py         BartlettArray class. Zero-endpoint triangular,
                          ~ -26 dB (contrast with triangular_array's nonzero
                          edges).
cosine_array.py           CosineArray class. Half-cosine taper, ~ -23 dB.
cosine_squared_array.py   CosineSquaredArray class. Raised cosine, ~ -31 dB.
hann_array.py             HannArray class. Identical to cosine-squared.
hamming_array.py          HammingArray class. 0.54/0.46 raised cosine, ~ -41 dB.
blackman_array.py         BlackmanArray class. Three-term cosine, ~ -58 dB.
kaiser_array.py           KaiserArray class. Parametric taper; shape parameter
                          beta trades beamwidth vs. sidelobe level.
woodward_lawson_array.py  WoodwardLawsonArray class. Shaped-beam synthesis
                          (flat-top, cosecant-squared) via composing beams.
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
                            bartlett_array,cosine_array,cosine_squared_array,
                            hann_array,hamming_array,blackman_array,kaiser_array,
                            dolph_tschebyscheff,taylor_array,villeneuve_array,
                            woodward_lawson_array,evaluate}
                           ...

Array Calculator

positional arguments:
  {uniform_array, binomial_array, triangular_array, bartlett_array,
   cosine_array, cosine_squared_array, hann_array, hamming_array,
   blackman_array, kaiser_array, dolph_tschebyscheff, taylor_array,
   villeneuve_array, woodward_lawson_array, evaluate}
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
  --norm {edge,center}  Normalize amplitudes to the edge or center element
                        (default center; dolph_tschebyscheff defaults to edge)
```

### CLI vs Scripting Comparison

ADD EXAMPLE HERE of the function calls in the CLI format, but just directly from python scripts
This is only typically used when integrating the calculators into other programs without losing the CLI functionality
but there's also some relevance because it opens up the calculator to scripting options. TODO, eventually. 

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


**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
shell = ArrayCalculator(['uniform_array', '-N', '10', '-f', '3e9'])
shell.main(shell.getArgs())
```
Results:

```bash

```
    

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Binomial array, verbose (also prints raw Pascal-triangle amplitudes) ---
    # python array_calculator.py binomial_array -N 6 -f 3e9 --verbose
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '-f', '3e9', '--verbose'])
    # shell.main(shell.getArgs())
    
```
Results:

```bash

```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Dolph-Tschebyscheff, -26 dB sidelobes, edge-normalized ---
    # python array_calculator.py dolph_tschebyscheff -N 10 -sll 26 --norm edge
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '10', '-sll', '26', '--norm', 'edge'])
    # shell.main(shell.getArgs())
```
Results:

```bash

```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Dolph-Tschebyscheff, center-normalized variant ---
    # python array_calculator.py dolph_tschebyscheff -N 12 -sll 30 --norm center
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '12', '-sll', '30', '--norm', 'center'])
    # shell.main(shell.getArgs())
```
Results:

```bash

```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Triangular taper, save a polar plot to a file (headless-safe) ---
    # python array_calculator.py triangular_array -N 16 --save tri.png --plot-style polar
    # shell = ArrayCalculator(['triangular_array', '-N', '16', '--save', 'tri.png', '--plot-style', 'polar'])
    # shell.main(shell.getArgs())

```
Results:

```bash

```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python

    # --- Steer the beam to 60 degrees (works on any array) ---
    # python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --scan 60 --verbose
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '8', '-sll', '25', '--scan', '60', '--verbose'])
    # shell.main(shell.getArgs())
```
Results:

```bash

```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Ordinary endfire: steer to the array axis (use closer spacing) ---
    # python array_calculator.py uniform_array -N 10 -d 0.25 --scan 0
    # shell = ArrayCalculator(['uniform_array', '-N', '10', '-d', '0.25', '--scan', '0'])
    # shell.main(shell.getArgs())
```
Results:

```bash

```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Export the radiation pattern to CSV (with steering) ---
    # python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --csv pattern.csv --scan 60
    # shell = ArrayCalculator(['dolph_tschebyscheff', '-N', '8', '-sll', '25', '--csv', 'pattern.csv', '--scan', '60'])
    # shell.main(shell.getArgs())

```
Results:

```bash


```

**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    # --- Return variables instead of printing, then read them back ---
    # python array_calculator.py binomial_array -N 6 --variable_return
    # shell = ArrayCalculator(['binomial_array', '-N', '6', '--variable_return'])
    # shell.main(shell.getArgs())
    # amps, hpbw_deg, directivity_db = shell.getCalcedParams()
    # print("done!")

```
Results:

```bash


```
**Uniform array (reference case), physical spacing from frequency**

CLI input: `python array_calculator.py uniform_array -N 10 -f 3e9`

Script:
```Python
    #   # evaluate returns a results dict, not a tuple:
    # shell = ArrayCalculator(['evaluate', '-g', './example_data/geometry.csv', '--variable_return'])
    # shell.main(shell.getArgs())
    # results = shell.getCalcedParams()        # results['peak_sidelobe_db'], etc.
    # print(results)

```
Results:

```bash


```
   
 

   







### Beam steering and endfire

Every array accepts `--scan`, the angle (in degrees) the main beam is steered
to. The default is `90` (broadside, beam perpendicular to the array axis).
Steering is applied as a uniform progressive phase across the elements,
`beta = -2*pi*(d/lambda)*cos(scan)`, so the same excitation amplitudes work at
any scan angle: the synthesis method (binomial, Dolph, the tapers, Taylor) sets
the *shape* of the beam, and `--scan` sets *where it points*.

```
# broadside (default)
python array_calculator.py uniform_array -N 10

# steered to 60 degrees
python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --scan 60

# ordinary endfire: steer all the way to the array axis
python array_calculator.py uniform_array -N 10 -d 0.25 --scan 0
```

**Endfire is just `--scan 0` (or `--scan 180`).** There is no separate "endfire
array" mode, because endfire is not a separate kind of array -- it is the
broadside synthesis steered to the array axis. Setting `--scan 0` (or `180`)
applies the ordinary-endfire progressive phase automatically, and any of the
amplitude distributions can be steered to endfire this way. For a usable endfire
beam you will generally want closer element spacing (`-d 0.25` or so); at
`d = lambda/2` the steered beam cannot fully reach the axis without its mirror
image at `180 - scan` merging into it, which distorts the pattern near endfire.

Two limitations are left to the user to account for, by design -- the calculator
does not detect or correct them:

1. **Directivity is the array-factor (isotropic-element) value and assumes a
   broadside-type pattern.** It is reported from the standard broadside
   expressions, which hold across scan angles for the isotropic-element model
   *until the beam approaches endfire*. Near `--scan 0` or `--scan 180` the
   broadside formula no longer applies (ordinary endfire is actually higher, on
   the order of `4 N d / lambda` for the uniform case), so the printed
   directivity at or near endfire should not be trusted. The radiation pattern,
   HPBW, and CSV/plot output remain correct at any scan; only the single
   directivity figure is model-limited.

2. **Real element patterns are not modeled.** The tool computes the array factor
   only, treating every element as isotropic. A real element radiates less
   toward endfire, which in practice reduces scanned directivity. If you need
   scan-loss or true scanned gain, fold in your element pattern separately.

The Hansen-Woodyard increased-directivity endfire condition is **not** provided
and is not reachable through `--scan` alone: it adds an extra phase increment
beyond the geometric endfire condition (`beta = -(2*pi*d/lambda + delta)`, with
`delta` on the order of `2.94/N`), deliberately over-steering past the axis to
trade higher sidelobes for greater directivity. That is a distinct synthesis
choice rather than a steering angle, so it is out of scope here; apply the extra
phase increment to the excitations yourself if you need it.

### Closed-form taper arrays

`uniform_array`, `binomial_array`, `triangular_array`, `bartlett_array`,
`cosine_array`, `cosine_squared_array`, `hann_array`, `hamming_array`, and
`blackman_array` take only the common arguments. For example:

```
python array_calculator.py triangular_array -N 16
python array_calculator.py bartlett_array -N 16
python array_calculator.py blackman_array -N 20 --plot --plot-style polar
python array_calculator.py hamming_array -N 16 --save hamming.png
```

### Dolph-Tschebyscheff Usage

Adds `-sll` (required) to the common arguments. The shared `--norm` flag still
applies, but unlike the other arrays it defaults to `edge` here, since
Dolph-Tschebyscheff coefficients are conventionally reported edge-normalized.

```
usage: array_calculator.py dolph_tschebyscheff [common args]
                                               -sll SIDELOBE_LEVEL
                                               [--norm {edge,center}]

  -sll, --sidelobe_level SIDELOBE_LEVEL
                        Desired sidelobe level in dB (e.g. 26, 30) (required)
  --norm {edge,center}  Normalize to edge or center element (default edge for
                        this subcommand; common arg, see above)
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

### Villeneuve Array Usage

Same arguments as Taylor (`-sll` required, `-nbar` optional). Villeneuve is the
discrete-element analog of Taylor: rather than sampling a continuous line-source
distribution, it places the pattern zeros directly for a finite N (stretched
Dolph-Tschebyscheff zeros inside `nbar`, uniform-array zeros beyond). For small
and moderate N this lands the near-in sidelobes on the design level more
accurately than sampled Taylor; the two converge as N grows.

```
usage: array_calculator.py villeneuve_array [common args]
                                            -sll SIDELOBE_LEVEL [-nbar NBAR]

  -sll, --sidelobe_level SIDELOBE_LEVEL
                        Desired (near-in) sidelobe level in dB (required)
  -nbar, --nbar NBAR    Number of near-in sidelobes held near the design level
                        (default 5)
```

### Kaiser Array Usage

Adds `-beta` to the common arguments. The single shape parameter `beta` trades
main-beam width against sidelobe level continuously: `beta = 0` is uniform
(highest sidelobes, narrowest beam) and larger `beta` lowers the sidelobes at the
cost of a wider beam. It is the closed-form taper closest to the prolate-
spheroidal optimum, and a useful bridge between the fixed tapers and the
parametric synthesis methods.

```
usage: array_calculator.py kaiser_array [common args] [-beta BETA]

  -beta, --beta BETA    Kaiser shape parameter (larger = lower sidelobes, wider
                        beam; default 6.0)
```

### Woodward-Lawson Usage

Woodward-Lawson is different from every other method here: instead of producing a
pencil beam with controlled sidelobes, it designs the array to *approximate a
desired pattern shape*. Two target shapes are built in -- `flat_top` (uniform
amplitude across a sector, zero outside) and `cosecant_squared` (the classic
radar ground-coverage pattern). The shaped sector is centered on `--scan` and its
half-width is set by `--sector`.

```
usage: array_calculator.py woodward_lawson_array [common args]
                                                 [--shape {flat_top,cosecant_squared}]
                                                 [--sector SECTOR] [--floor FLOOR]

  --shape {flat_top,cosecant_squared}
                        Target pattern shape (default flat_top)
  --sector SECTOR       Half-width of the shaped sector in degrees (default 30)
  --floor FLOOR         Desired linear amplitude outside the sector (default 0)
```

Because the excitations are generally complex, `--verbose` also prints the
per-element excitation phase. The in-sector ripple is an inherent (Gibbs-like)
feature of the composing-beam method and decreases as N grows.

### Evaluating an arbitrary geometry

The `evaluate` subcommand is not a synthesis method -- it *scores* an array you
already have. Given a geometry CSV with element positions (and optionally
per-element amplitudes and phases), it runs the full analysis pipeline (pattern
sweep, main-beam angle, HPBW, peak sidelobe, directivity-by-integral) and can
export or plot the result like any other array. This is the intended boundary
with an external optimizer: the optimizer proposes a geometry, the evaluator
scores it. No optimizer lives in the calculator.

```
usage: array_calculator.py evaluate [common args] -g GEOMETRY_CSV

  -g, --geometry GEOMETRY_CSV
                        Geometry CSV with columns:
                        position_lambda[, amplitude][, phase_deg]
                        (amplitude defaults to 1.0, phase_deg to 0.0)
```

The geometry CSV looks like:

```
position_lambda,amplitude,phase_deg
0.0,1.0,0
0.5,1.0,0
1.0,1.0,0
...
```

Programmatically (for an optimizer's fitness call), use the method directly:

```python
from array_evaluator import ArrayEvaluator
ev = ArrayEvaluator()
result = ev.evaluate(positions, amplitudes=None, phases=None)
# result is a dict: peak_theta_deg, hpbw_deg, peak_sidelobe_db,
# directivity, directivity_db, plus the raw pattern arrays.
```

Directivity here is computed from the pattern integral (valid for non-uniform
spacing), not the `2 N d / lambda` closed form, which only applies to equally
spaced arrays.

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
[*] Directivity = 9.19
[*] Directivity = 9.63 dB
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

**Villeneuve n-bar array, 20 elements, -30 dB near-in sidelobes, nbar = 6:**
```
python array_calculator.py villeneuve_array -N 20 -sll 30 -nbar 6

[*] N = 20
[*] Sidelobe level = 30.00 dB
[*] nbar = 6
[*] Amplitudes (normalized) = [0.266, 0.301, 0.378, 0.490, 0.612, 0.725, 0.826, 0.910, 0.970, 1.000, 1.000, 0.970, 0.910, 0.826, 0.725, 0.612, 0.490, 0.378, 0.301, 0.266]
[*] Element spacing d = 0.50 lambda
[*] HPBW = 6.43 deg
[*] Directivity = 17.15
[*] Directivity = 12.34 dB
```

Villeneuve is the discrete analog of Taylor; for this N the two distributions are
nearly identical (compare the amplitudes with the Taylor example above), but
Villeneuve places the near-in sidelobes on the design level more precisely as N
shrinks.

**Bartlett taper, 16 elements:**
```
python array_calculator.py bartlett_array -N 16

[*] N = 16
[*] Amplitudes (normalized) = [0.000, 0.143, 0.286, 0.429, 0.571, 0.714, 0.857, 1.000, 1.000, 0.857, 0.714, 0.571, 0.429, 0.286, 0.143, 0.000]
[*] Element spacing d = 0.50 lambda
[*] HPBW = 9.82 deg
[*] Directivity = 11.20
[*] Directivity = 10.49 dB
```

Bartlett is the zero-endpoint triangular taper: note the edge elements go to
exactly 0, unlike `triangular_array`, which keeps them nonzero.

**Kaiser taper, 16 elements, beta = 8:**
```
python array_calculator.py kaiser_array -N 16 -beta 8

[*] N = 16
[*] beta = 8
[*] Amplitudes (normalized) = [0.002, 0.027, 0.096, 0.231, 0.430, 0.663, 0.874, 1.000, 1.000, 0.874, 0.663, 0.430, 0.231, 0.096, 0.027, 0.002]
[*] Element spacing d = 0.50 lambda
[*] HPBW = 12.13 deg
[*] Directivity = 9.01
[*] Directivity = 9.55 dB
```

Raising `beta` lowers the sidelobes and widens the beam; lowering it toward 0
approaches the uniform array.

**Woodward-Lawson flat-top beam, 20 elements, +/-30 deg sector:**
```
python array_calculator.py woodward_lawson_array -N 20 --shape flat_top --sector 30

[*] N = 20
[*] Target shape = flat_top
[*] Amplitudes (normalized) = [0.079, 0.081, 0.085, 0.092, 0.103, 0.121, 0.150, 0.205, 0.336, 1.000, 1.000, 0.336, 0.205, 0.150, 0.121, 0.103, 0.092, 0.085, 0.081, 0.079]
[*] Element spacing d = 0.50 lambda
[*] Directivity = 8.24
[*] Directivity = 9.16 dB
```

Woodward-Lawson approximates a *shape* rather than a pencil beam; here the
pattern is held roughly flat across a 60-degree sector and suppressed outside it.
Use `--shape cosecant_squared` for the radar ground-coverage pattern, and
`--verbose` to also print the per-element excitation phase.

**Evaluate an arbitrary geometry from a CSV:**
```
python array_calculator.py evaluate -g geometry.csv

[*] N = 8
[*] Positions (lambda) = [0.000, 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500]
[*] Amplitudes = [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000]
[*] Main beam at = 90.00 deg
[*] HPBW = 12.81 deg
[*] Peak sidelobe = -12.80 dB
[*] Directivity = 8
[*] Directivity = 9.03 dB
```

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

The test suite uses `pytest` and is split one file per array
(`test_uniform_array.py`, `test_binomial_array.py`, `test_dolph_tschebyscheff.py`,
and so on), with shared fixtures in `test_helpers.py` and cross-cutting checks
(array factor / steering, CSV export, CLI dispatch, edge cases) in
`test_array_calculator.py`. It checks the synthesis math against the worked
examples and verifies the computed patterns physically (e.g. that
a Dolph-Tschebyscheff design's sidelobes actually land at the requested level,
and that a binomial array has none).

`test_helpers.py` defines only the shared helpers (`make_args`,
`peak_sidelobe_db`) and is imported by the other files; despite the `test_`
prefix it collects no tests of its own.

Run the full suite with:

```
pytest -v
```

Or a single array's tests directly:

```
pytest -v tests/test_triangular_array.py
```

Coverage includes:
* Binomial amplitudes (Pascal's triangle), HPBW, and directivity
* Dolph-Tschebyscheff z0, R, and coefficient synthesis for even and odd N
* Villeneuve: peak sidelobe lands on the design level across N/sll/nbar
* Edge- vs. center-normalization: the `--norm` toggle changes only the printed
  scale, never the ordering or the radiation pattern (sidelobe level invariant)
* Dolph-Tschebyscheff directivity tracks the directivity integrated from the
  actual pattern (guards against the broadside-length scaling bug)
* Closed-form tapers: symmetry, center peak, and the characteristic peak
  sidelobe level of each (triangular / Bartlett ~ -26 dB, cosine ~ -23 dB,
  cosine-squared / Hann ~ -31 dB, Hamming ~ -41 dB, Blackman ~ -58 dB)
* Kaiser: I0 series vs. reference values; sidelobes fall monotonically with beta
* The Hann / cosine-squared equivalence (identical amplitudes across N)
* Woodward-Lawson: the realized pattern approximates the requested shape
  (flat-top holds its sector and drops outside; cosecant-squared is broad)
* Evaluator round-trip: a uniform geometry scored by the evaluator matches the
  dedicated UniformArray; phase steers the beam; non-uniform spacing runs cleanly
* Physical pattern checks: target sidelobe level, equiripple sidelobes,
  no-sidelobe binomial behavior, and beam steering
* CSV export format (pattern export and geometry import)
* CLI dispatch, defaults, `--variable_return`, and argument validation

Running `pytest` requires `pytest` in addition to the project's `numpy` and
`pint` dependencies.


## Architecture and Flow
 
Addition for future work and integration with other software. This project is primarily a calculator, but the evaluator is useful for recreating and comparing results with the same suite. 

There's two ways to input data, and a shared engine for the output. The **synthesis** path *generates* an
excitation from a design goal; the **evaluator** path *accepts* an arbitrary geometry from outside (e.g. an optimizer). Both feed the same `ArrayCommon`
pipeline that turns excitations into a pattern, figures of merit, and output.
 
```
                              ┌───────────────────────────┐
                              │         USER / CLI         │
                              │   array_calculator.py      │
                              │   (argparse dispatch)      │
                              └─────────────┬─────────────┘
                                            │
                ┌───────────────────────────┴───────────────────────────┐
                │                                                         │
                ▼                                                         ▼
   ════════ SYNTHESIS PATH ════════                          ════════ EVALUATE PATH ════════
   "design an array from a goal"                             "score a geometry I already have"
                │                                                         │
                ▼                                                         ▼
   ┌─────────────────────────────┐                          ┌─────────────────────────────┐
   │  one synthesis subcommand   │                          │   evaluate  -g geometry.csv │
   │ ─────────────────────────── │                          │ ─────────────────────────── │
   │ uniform   binomial          │                          │  ArrayEvaluator             │
   │ dolph_t.  taylor            │                          │  read_geometry_csv()        │
   │ villeneuve                  │                          │   position_lambda           │
   │ triangular bartlett         │                          │   amplitude   (opt)         │
   │ cosine  cosine_squared      │                          │   phase_deg   (opt)         │
   │ hann  hamming  blackman     │                          └──────────────┬──────────────┘
   │ kaiser                      │                                         │
   │ woodward_lawson             │            (optimizer can call)         │
   └──────────────┬──────────────┘            evaluate(positions,          │
                  │                            amplitudes, phases) ───────►│
                  ▼                                                         │
   ┌─────────────────────────────┐                                         │
   │ amplitudes()  [+ phases]    │   produces excitation                   │
   │ (per-method synthesis math) │   (amps, and phases for                 │
   └──────────────┬──────────────┘    woodward_lawson)                     │
                  │                                                         │
                  │   excitation = positions + amplitudes + phases          │
                  └─────────────────────────┬───────────────────────────────┘
                                            │
                                            ▼
                       ┌─────────────────────────────────────┐
                       │            ArrayCommon              │
                       │   (shared engine — one source of    │
                       │        truth for the physics)       │
                       │ ─────────────────────────────────── │
                       │  array_factor(amps, theta,          │
                       │               positions, phases)    │
                       │        AF(θ) = Σ aₙ e^{j(2π xₙcosθ  │
                       │                        + φₙ)}        │
                       │                 │                   │
                       │                 ▼                   │
                       │  pattern_sweep()  → |AF| over θ     │
                       │                 │                   │
                       │      ┌──────────┼──────────┐        │
                       │      ▼          ▼          ▼        │
                       │   HPBW    peak sidelobe  directivity│
                       │           (−3 dB)   (local maxima) (∫ or 2Nd/λ)│
                       └─────────────────────┬───────────────┘
                                            │
                       ┌─────────────────────┼───────────────────┐
                       ▼                     ▼                   ▼
              ┌────────────────┐   ┌──────────────────┐  ┌────────────────┐
              │  print to CLI  │   │  export CSV      │  │  plot pattern  │
              │  [*] N, amps,  │   │  theta_deg,      │  │  polar / rect  │
              │  HPBW, D, ...  │   │  AF_linear,AF_dB │  │  (matplotlib)  │
              └────────────────┘   └──────────────────┘  └────────────────┘
                       │
                       │  --variable_return
                       ▼
              ┌────────────────────────────┐
              │ getCalcedParams()          │
              │  synthesis → tuple         │
              │  evaluate  → results dict ─┼──►  external optimizer
              └────────────────────────────┘     (proposes next geometry,
                                                   loops back to EVALUATE PATH)
```


## Development and future work

* **Non-uniform / thinned / aperiodic synthesis** is intentionally *not* built in.
  Those are optimization problems and belong in the companion optimizer project;
  the `evaluate` subcommand is the hook for it (the optimizer proposes element
  positions, this calculator scores them). Statistical thinning and minimum-
  redundancy placement can be layered on top of the evaluator later without
  changing it.
* **Hansen-Woodyard** increased-directivity endfire is not provided as a mode
  (see [Beam steering and endfire](#beam-steering-and-endfire)); it can still be
  *evaluated* by passing the appropriate per-element phases to the evaluator.
* **Directivity near endfire** is reported from the broadside formula and is not
  valid as the scan angle approaches the array axis (see the steering section).
* **Woodward-Lawson defaults** (sector, floor, the in-sector ripple vs. N
  tradeoff) are reasonable but not yet tuned against worked teaching examples.
* The staged `windowed_array.py` consolidation is not yet wired into the CLI.

## References

[1]: C. A. Balanis, Antenna Theory: Analysis and Design, 4th ed. Hoboken, NJ, USA: Wiley, 2016.

[2]: V. Anand, "Antenna array factor calculations and plot," Vinoth.org. [Online]. Available: https://www.vinoth.org/rf-theory/antenna-array-factor-calculations-plot. Accessed: Jun. 4, 2026.

[3]: P. Delos, B. Broughton, and J. Kraft, "Phased array antenna patterns—Part 1: Linear array beam characteristics and array factor," Analog Dialogue, Analog Devices, 2020. [Online]. Available: https://www.analog.com/en/resources/analog-dialogue/articles/phased-array-antenna-patterns-part1.html. Accessed: Jun. 4, 2026.

[4]: "Uniform linear arrays," Antenna-Theory.com. [Online]. Available: https://www.antenna-theory.com/arrays/weights/uniform.php. Accessed: Jun. 4, 2026.

[5]: S. W. Ellingson, "100: Angle-of-arrival estimation (Bartlett)," YouTube, 2022. [Online video]. Available: https://www.youtube.com/watch?v=GC4ETN3Lhzk. Accessed: Jun. 4, 2026.

[6]: I. Berrios, "Introduction to beamforming, Part 2: The Bartlett beamformer," Medium, Sep. 30, 2024. [Online]. Available: https://medium.com/@itberrios6/introduction-to-beamforming-part-2-68db43c073b6. Accessed: Jun. 4, 2026. (and some of the other parts, but mostly this one)

[7]: "Antenna design and analysis, Session 32," Centurion University of Technology and Management, Courseware. [Online]. Available: https://courseware.cutm.ac.in/wp-content/uploads/2020/05/Antenna-Design-Analysis-Session-32-1.pdf. Accessed: Jun. 4, 2026.

[8]: Engineering Funda, "Binomial array (basics, pattern multiplication, Pascal's triangle & parameters) explained," YouTube. [Online video]. Available: https://youtu.be/tOkUtxvjvXA?list=PLgwJf8NK-2e6xvkHGQDZhqvHIW2BV0EpN. Accessed: Jun. 4, 2026.

[9]: S. V. Hum, "Antenna arrays II," ECE422 course notes, Univ. of Toronto, Toronto, ON, Canada. [Online]. Available: https://www.waves.utoronto.ca/prof/svhum/ece422/notes/15-arrays2.pdf. Accessed: Jun. 4, 2026.

[10]: "Array pattern synthesis," MathWorks Phased Array System Toolbox Documentation, The MathWorks, Inc. [Online]. Available: https://www.mathworks.com/help/phased/ug/array-pattern-synthesis.html. Accessed: Jun. 4, 2026.

[11]: H. Mistialustina, Chairunnisa, and A. Munir, "Evaluation of Kaiser function-based linear array performance in suppressing SLL and its experimental approach," IEEE Access, vol. 12, pp. 94712–94732, 2024, doi: 10.1109/ACCESS.2024.3424237.

[12]: A. Kumar, "Lecture 2 | Woodward-Lawson method | Array synthesis | Antenna and wave propagation," YouTube, Dr. Ashok Kumar. [Online video]. Available: https://www.youtube.com/watch?v=7uHRTEsYqe0. Accessed: Jun. 4, 2026.

[13]: H. Mistialustina, Chairunnisa, and A. Munir, "Analytical approach to parameter determination in Kaiser function for power-weighted antenna array design," J. ICT Res. Appl., vol. 17, no. 1, pp. 113–127, 2023, doi: 10.5614/itbj.ict.res.appl.2023.17.1.8.