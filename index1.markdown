---
layout: default
title: Effective Mass Calculator
---

## {{ page.title }}

[![Build Status](https://travis-ci.org/alexandr-fonari/emc.png)](https://travis-ci.org/alexandr-fonari/emc.png)

### Theory

Effective mass calculator (**EMC**) implements calculation of the effective masses at the bands extrema using finite difference method. Effective mass (m*) is defined as:

![Eq. 1](/emc/eqs/01.svg)

where *x, y, z* are the directions in the reciprocal Cartesian space (2π/A), *En(k)* is the dispersion relation for the *n*-th electronic band. The explicit form of the right-side symmetric tensor from the above eq. is:

![Eq. 1](/emc/eqs/02.svg)

where the derivatives can be evaluated numerically using the finite difference method. Eigenvalues of the above matrix are inverses of the effective masses, eigenvectors are the directions of the principal effective mass components.

For the theory behind the code and validation of the code against known data see [the paper](https://github.com/alexandr-fonari/emc/blob/master/Paper-03-18-2013.pdf?raw=true). Let us know if you find any bugs or mistakes, thanks!

#### Notes
 1. Atomic units (a.u.) are used throughout the code: hbar = 1, energy is in Hartree, distance is in Bohr, mass is in the electron mass at rest (m0).
 1. For the top of the VB (valence band) eigenvalues are negative, for the bottom of the CB (conduction band) eigenvalues are positive (as results from the basic calculus).
 1. In some cases, not all eigenvalues have the same sign, meaning that the chosen k-point is not a global minimum (maximum).
 1. Effective masses can be highly anisotropic.

### Input file structure

```bash
0.000 0.000 0.000                       ! K-POINT in the reciprocal crystal coord. (3 floats)
0.01                                    ! step size in 1/Bohr units (1 float)
81                                      ! band number, (1 integer)
V                                       ! program identifier (1 char)
6.291999817  0.000000000  0.000000000   ! direct lattice vectors (3 floats)
0.755765092  7.652872670  0.000000000   ! direct lattice vectors (3 floats)
0.462692761  3.245907103 14.032346772   ! direct lattice vectors (3 floats)
```

 - **band number**. If *CRYSTAL* is employed, band number should be set to **1**. Helper script `cry-getE.pl` reads-in the desired band number (see below). For *VASP* valence band number can be obtained as a half of the `NELECT` variable from the *OUTCAR* file (for non spin-polarized calculations).
 - **program identifier**: `C` for *CRYSTAL* or `V` for *VASP*. (TODO: Quantum Espresso)
 - **direct lattice components** in *CRYSTAL* can be found under: `DIRECT LATTICE VECTORS COMPON. (A.U.)`, in *VASP* under: `direct lattice vectors`.

### 3. Usage
#### 3.1 CRYSTAL
1. Run SCF.
1. Create a directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```.
1. Create ```inp``` file with the desired characteristics in the newly created directory.
1. Run ```emc_gen.x``` to obtain ```KPOINTS``` file.
1. Run ```cry-getE.pl -f ../input.f9 -b 131``` to obtain ```EIGENVAL``` file.  
 **cry-getE.pl** uses the following parameters:
   * SCF output filename (.f9) in ```-f``` flag;
   * band number in ```-b``` flag.
   * Note that ```runprop09``` needs to be in the current ```$PATH```, otherwise script will quit.
1. Run ```emc_calc.x``` to obtain effective masses and directions. Look for ```emc_calc.log``` file.

#### 3.2 VASP
1. Run SCF (e.g. ```ICHARG=2``` in INCAR).
1. Create a directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```.
1. Create ```inp``` file with the desired characteristics in the newly created directory.
1. Run ```emc_gen.x``` to obtain ```KPOINTS``` file.
1. Run a non-SCF calculation (```ICHARG=11``` in INCAR). Don't forget to copy ```CHGCAR``` file from SCF folder.
1. Run ```emc_calc.x``` to obtain effective masses and directions. Look for ```emc_calc.log``` file.

### 4. Acknowledgments and references
1. Finite-difference method on three-point stencil is outlined here: K. Doi, *et al.*, *J. Appl. Phys.*, **98**, 113709 (2005): [10.1063/1.2138381](http://dx.doi.org/10.1063/1.2138381).

### 5. Test cases
In progress...