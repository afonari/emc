# Effective Mass Calculator for Semiconductors

Effective mass calculator (**EMC**) implements calculation of the effective masses at the bands extrema using finite difference method (**not** the band fitting method). There are two versions of the program: written in FORTRAN and Python. Currently *CRYSTAL*, *VASP*, *CASTEP* are supported, *Quantum Espresso* is coming!

[![Build Status](https://travis-ci.org/afonari/emc.svg?branch=master)](https://travis-ci.org/afonari/emc) [![Coverage Status](https://coveralls.io/repos/github/afonari/emc/badge.png?branch=master)](https://coveralls.io/github/afonari/emc?branch=master)

## Theory

Effective mass (m*) is defined as:

![Eq. 1](https://rawgithub.com/afonari/emc/master/images/01.svg)

where *x, y, z* are the directions in the reciprocal Cartesian space (2π/A), *En(k)* is the dispersion relation for the *n*-th electronic band. The explicit form of the right-side symmetric tensor from the above eq. is:

![Eq. 1](https://rawgithub.com/afonari/emc/master/images/02.svg)

where the derivatives can be evaluated numerically using the finite difference method. Eigenvalues of the above matrix are inverses of the effective masses, eigenvectors are the directions of the principal effective mass components.

For the theory behind the code and validation of the code against known data see [the paper](https://github.com/alexandr-fonari/emc/blob/master/Paper-03-18-2013.pdf?raw=true). Let us know if you find any bugs or mistakes, thanks!

#### Notes
 1. Atomic units (a.u.) are used throughout the code: ħ = 1, energy is in Hartree, effective mass is in the electron mass at rest (m0), distance is in Bohr
 1. For the top of the VB (valence band) eigenvalues are negative, for the bottom of the CB (conduction band) eigenvalues are positive
 1. In some cases, not all eigenvalues have the same sign, meaning that the chosen k-point is not a global minimum (maximum)
 1. Effective masses can be highly anisotropic (see [tests](#toc_11) section)

## Installation

### Python version

`emc.py` is a Python script, that depends only on the Python Standard Library. Code is being tested with Python v 2.7. 

To install:
 - check that *emc.py* has executable flag using `ls -la`, if it doesn't, do `chmod +x ./emc.py`
 - check that *emc.py* is in your path `$PATH` (to print the `$PATH` variable do `echo $PATH`)
 - enjoy the results!

### Fortran version (considered deprecated)

FORTRAN version is considered deprecated: Download and unpack the current version: [**1.50**](https://github.com/alexandr-fonari/emc/releases/download/1.50/emc-1.50.tar.gz).

Fortran version depends on LAPACK and BLAS libraries.

To install:

 - edit Makefile from the *fortran* folder for you needs
 - run `make`
 - check that *emc_gen* and *emc_calc* are in your path `$PATH` (to print the `$PATH` variable do `echo $PATH`)
 - enjoy the results!

NOTE: `emc.py` uses STENCIL=3 by default. It is possible to use STENCIL=5. For that purpose change the third line in `emc.py` file from `STENCIL=3` to `STENCIL=5`.

## Input file structure

```bash
0.000 0.000 0.000                       ! K-POINT in the reciprocal crystal coord. (3 floats)
0.01                                    ! step size in 1/Bohr units (1 float)
81                                      ! band number, (1 integer)
V                                       ! program identifier (1 char)
6.291999817  0.000000000  0.000000000   ! direct lattice vectors (3 floats)
0.755765092  7.652872670  0.000000000   ! direct lattice vectors (3 floats)
0.462692761  3.245907103 14.032346772   ! direct lattice vectors (3 floats)
```

 - **band number**. If *CRYSTAL* is employed, band number should be set to **1**. Helper script `cry-getE.pl` reads-in the desired band number ([see below](#toc_8)). For *VASP* valence band number can be obtained as a half of the `NELECT` variable from the *OUTCAR* file (for non spin-polarized calculations).
 - **program identifier**: `C` for *CRYSTAL* or `V` for *VASP*. (TODO: Quantum Espresso)
 - **direct lattice components** in *CRYSTAL* can be found in the SCF output under: `DIRECT LATTICE VECTORS COMPON. (A.U.)`, in *VASP* in the *OUTCAR* under: `direct lattice vectors`.

## Usage

1. Run SCF.
1. Generate k-point grid (in *KPOINT* file):
    - Python version: `emc.py input_file`
    - FORTRAN version: `emc_gen`. Note: FORTRAN version requires input file to be named `inp`
1. Run non-self consistent calculation using obtained k-point grid:
    - in the case of CRYSTAL: run helper script cry-getE.pl (see below)
    - in the case of VASP: set `ICHARG=11` in the *INCAR*. Don't forget to copy *CHGCAR* file from the converged SCF run
1. Calculate effective masses and principal directions using *EIGENVAL* file generated in previous step:
    - Python version: `emc.py input_file EIGENVAL_file`
    - FORTRAN version: `emc_calc`. Current folder should contain both *inp* and *EIGENVAL* files

### Helper script for CRYSTAL: cry-getE.pl

In case of *CRYSTAL*, *cry-getE.pl* script should be used in order to obtain file with the energies on the grid. The script takes two k-points at a time and runs band structure calculations (using *runprop* script from the *CRYSTAL* package).

*cry-getE.pl* has the following command line options:

 - `-f` SCF output filename (.f98)
 - `-b` band number

Example: `cry-getE.pl -f ../input.f98 -b 131`

Note that ```runprop``` needs to be in the current ```$PATH```, otherwise script will quit.

## Usage with CASTEP
*Contribution by Genadi Naydenov*

1. Create an input file (let's say emc-input). Below is an example of emc-input:

```bash
0.000 0.000 0.000                       ! K-POINT in the reciprocal crystal coord. (3 floats)
0.01                                    ! step size in 1/Bohr units (1 float)
17                                      ! band number, (1 integer)
P                                       ! program identifier (1 char) - P is the CASTEP identifier
6.291999817  0.000000000  0.000000000   ! direct lattice vectors in Angs (3 floats)
0.755765092  7.652872670  0.000000000   ! direct lattice vectors in Angs (3 floats)
0.462692761  3.245907103 14.032346772   ! direct lattice vectors in Angs (3 floats)
```
2. Run `./emc.py emc-input`. This will generate a `KPOINTS` file, which contains a list of k-points.
3. Copy the fractional coordinates of the k-points from KPOINTS and paste them into a kpoints_list block in <seedname>.cell. Example of the block:

```bash
%BLOCK SPECTRAL_KPOINTS_LIST

*paste k-points list here*

%ENDBLOCK SPECTRAL_KPOINTS_LIST
```

4. Run a spectral calculation in CASTEP.
5. Once the spectral calculation is done, use `./emc.py emc-input <seedname>.bands` to obtain the effective mass for the desired band. You can change the band number and repeat step 5 to calculate the effective mass for a different band.

## Running tests

To run tests, in the distribution directory run:  
`python -m unittest discover`

Tests are located in the *test* folder.

## Authors

Alexandr Fonari and Christopher Sutton  
If you have any questions or suggestions don't hesitate to contact us at: alexandr[dot]fonari[nospam]gmail.com or csutton[nospam]gatech.edu. You can also [submit as issue](https://github.com/alexandr-fonari/emc/issues/new).

## License: MIT

Copyright (c) 2012, Alexandr Fonari, Christopher Sutton  
Cite as: "Effective Mass Calculator, A. Fonari, C. Sutton, (2012)."

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
