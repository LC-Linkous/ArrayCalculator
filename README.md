# ArrayCalculator
in progress AntennaCAT compatible array calculator.

## NOT STABLE, NOT COMPLETE


This project is a CLI-based linear antenna array synthesis (sort of, full functionality in progress) tool. 
It is written as an AntennaCAT [ADD LINK] compatible tool, similar to [ADD LINK] [AntennaCalculator]().
This tool designs the excitation coefficients for a multi-element array given a design
goal, then reports the beam characteristics and (optionally) exports the radiation pattern.



The calculator features the following synthesis methods:
* Binomial array (no sidelobes, amplitudes from Pascal's triangle)
* Dolph-Tschebyscheff array (narrowest beam for a specified sidelobe level)
* 2-3 more, including Taylor to be added from class notes 


Supported outputs:
* Excitation coefficients (edge- or center-normalized), HPBW, directivity
* Radiation-pattern export to CSV
* Radiation-pattern plot (optional, requires matplotlib)



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



## File Structure

```python

ArrayCalculator/
├── .python-version
├── pyproject.toml
├── README.md
├── LICENSE
├── .gitignore
├── src/
│   ├── __init__.py
│   ├── array_calculator.py
│   ├── array_common.py
│   ├── binomial_array.py
│   └── dolph_tschebyscheff.py
└── tests/
    ├── __init__.py
    ├── test_array_calculator.py
    └── ADD.py
```

```
array_calculator.py       Main entry point. Builds the CLI and dispatches
                          to the synthesis classes (ArrayCalculator class).
array_common.py           Shared base class (ArrayCommon): printing helpers,
                          frequency/spacing reporting, array-factor
                          computation, CSV export, and the plot hook.
binomial_array.py         BinomialArray class. Pascal's-triangle amplitudes,
                          HPBW, and directivity.
dolph_tschebyscheff.py    DolphTschebyscheff class. Coefficient synthesis for
                          a specified sidelobe level, with edge/center
                          normalization.
```

Both array classes inherit from `ArrayCommon`, so the two methods only supply
their excitation amplitudes; the shared pattern and output process lives in
one place. 

## Usage

```
usage: array_calculator.py [--help] [--version] {binomial_array,dolph_tschebyscheff} ...

Array Calculator

positional arguments:
  {binomial_array,dolph_tschebyscheff}
                        sub-command help

optional arguments:
  --help                Show this help message and exit
  --version             show program's version number and exit
```

### Binomial Array Usage
```
usage: array_calculator.py binomial_array [--help] [--verbose] -N ELEMENTS [-f FREQUENCY]
                                          [-d SPACING] [--scan SCAN]
                                          [-u {meter,centimeter,millimeter,inch}]
                                          [--csv CSV] [--plot] [--variable_return]

optional arguments:
  --help                Show this help message and exit
  --verbose
  -N ELEMENTS, --elements ELEMENTS
                        Number of array elements
  -f FREQUENCY, --frequency FREQUENCY
                        Frequency in Hz (optional; gives physical spacing)
  -d SPACING, --spacing SPACING
                        Element spacing as a fraction of wavelength (default 0.5)
  --scan SCAN           Beam steering angle in degrees (default 90 = broadside)
  -u {meter,centimeter,millimeter,inch}, --unit {meter,centimeter,millimeter,inch}
                        Unit of measurement
  --csv CSV             Export radiation pattern to CSV file
  --plot                Plot radiation pattern (requires matplotlib)
  --variable_return     Return variables instead of printing
```

### Dolph-Tschebyscheff Usage
```
usage: array_calculator.py dolph_tschebyscheff [--help] [--verbose] -N ELEMENTS -sll SIDELOBE_LEVEL
                                               [--norm {edge,center}] [-f FREQUENCY] [-d SPACING]
                                               [--scan SCAN] [-u {meter,centimeter,millimeter,inch}]
                                               [--csv CSV] [--plot] [--variable_return]

optional arguments:
  --help                Show this help message and exit
  --verbose
  -N ELEMENTS, --elements ELEMENTS
                        Number of array elements
  -sll SIDELOBE_LEVEL, --sidelobe_level SIDELOBE_LEVEL
                        Desired sidelobe level in dB (e.g. 26, 30)
  --norm {edge,center}  Normalize to edge or center element (default edge)
  -f FREQUENCY, --frequency FREQUENCY
                        Frequency in Hz (optional; gives physical spacing)
  -d SPACING, --spacing SPACING
                        Element spacing as a fraction of wavelength (default 0.5)
  --scan SCAN           Beam steering angle in degrees (default 90 = broadside)
  -u {meter,centimeter,millimeter,inch}, --unit {meter,centimeter,millimeter,inch}
                        Unit of measurement
  --csv CSV             Export radiation pattern to CSV file
  --plot                Plot radiation pattern (requires matplotlib)
  --variable_return     Return variables instead of printing
```

## Examples

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



**Export the pattern to CSV and steer the beam to 60 degrees:**
```
python array_calculator.py dolph_tschebyscheff -N 8 -sll 25 --csv pattern.csv --scan 60
```

The CSV contains three columns: `theta_deg`, `AF_linear`, and `AF_dB`.


[INSERT DATA EXAMPLE]



**Matplotlib visualization:**
```
python array_calculator.py ARGS
```

[INSERT IMGS]


This format is being finalized.


## Tests

The test suite is in `test_array_calculator.py` and uses `pytest`. It checks the
synthesis math against the worked examples in [1] and verifies the computed
patterns physically (e.g. that a Dolph-Tschebyscheff design's sidelobes actually
land at the requested level, and that a binomial array has none).

Run the tests with:

```
pytest -v
```

Or directly:

```
python test_array_calculator.py
```

Coverage includes:
* Binomial amplitudes (Pascal's triangle), HPBW, and directivity
* Dolph-Tschebyscheff z0, R, and coefficient synthesis for even and odd N
* Edge- vs. center-normalization consistency
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