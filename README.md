## Effective Mass Calculator

### 1. Theory

```
Note: Atomic units (a.u.) will be used, hbar = 1, energy is in Hartrees, distance is in Bohrs, mass is in electron mass (m0).
```

The average velocity of an electron for a certain value of *n* and *k* is:  
![Mean velocity](https://raw.github.com/alexandr-fonari/emc/master/pics/p_v.png)  
thus acceleration can be written as:  
![Acceleration](https://raw.github.com/alexandr-fonari/emc/master/pics/p_a.png)

In the semiclassical picture, equation of motion of an electron in the presence of the electric (**E**) and magnetic (**H**) fields is:  
![Equation of motion](https://raw.github.com/alexandr-fonari/emc/master/pics/p_e_m.png)  
in the case of: ```H = 0```:  
![Equation of motion](https://raw.github.com/alexandr-fonari/emc/master/pics/p_f.png)  
getting back to acceleration:  
![Acceleration](https://raw.github.com/alexandr-fonari/emc/master/pics/p_a2.png)  
from here, inverse of the effective mass tensor (9 comoponents) can be written as:  
![Inverse EM tensor](https://raw.github.com/alexandr-fonari/emc/master/pics/p_1o_m.png)  
Note that this tensor is symmetric (can be diagonalized with [DSYEV](http://netlib.org/lapack/double/dsyev.f)):  
![Tensor](https://raw.github.com/alexandr-fonari/emc/master/pics/p_tensor.gif)  
At a saddle point (e.g. band maximum/minimum) components of the effective mass are inverse of eigenvalues of the tensor:


Energy gradient tensor (symmetric) is defined in 3D as follows:  
![Energy Gradient Tensor](https://raw.github.com/alexandr-fonari/emc/master/pics/p_et.png)  
where ```x*, y*, z*``` are cartesian reciprocal directions (see ```EMCcoords.pl``` for reference).

Second order derivatives are estimated with ```O(h^4)``` error:  
![2nd Derivative](https://raw.github.com/alexandr-fonari/emc/master/pics/p_2ndd.png)

Mixed second derivatives are estimated also with ```O(h^4)``` error:  
![Mixed 2nd Derivative](http://www.holoborodko.com/pavel/wp-content/ql-cache/quicklatex.com-ead43440eddb0f8db2cc36a1df79c547_l3.svg)

#### 2. Required input files
**1.** ```OUTCAR``` can be either OUTCAR from an SCF calculation from VASP or output from an SCF calculation from CRYSTAL renamed to OUTCAR  
**2.** ```inp``` has the following form:  
```
0.000 0.000 0.000   ! K-POINT in reciprocal cartesian system, set to Gamma: 3 floats
0.001               ! dk step: 1 float
81                  ! band number, can be VB or CB or whatever you want: 1 integer
V                   ! program, currently support V for VASP and C for crystal: 1 char
```

#### 3. How to run with CRYSTAL
1. Run SCF.
1. Make directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```*.
1. Copy ```KPOINTS``` file generated with ```EMCg.x``` into it.
1. Run from freshly created directory (```emH-00-50-00-d01```) as ```cry-getE.pl -i ../input.out -f ../input.f9 -b 131```.
1. Copy ```EIGENVAL``` generated file to the directory with ```inp``` file.
1. Run ```EMCc.x``` in directory with ```inp``` and ```EIGENVAL``` files.
1. Determine real space directions from eigevectors using ```EMCcoords.pl``` script.  
*. In CRYSTAL, units are 1/Bohr.

#### 4. How to run with VASP
1. Run SCF (e.g. ```ICHARG=2``` in INCAR).
1. Make directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```*.
1. Copy ```KPOINTS``` file generated with ```EMCg.x``` into it.
1. Run non-SCF calculation (```ICHARG=11``` in INCAR) from freshly created directory (```emH-00-50-00-d01```)
1. Copy resulting ```EIGENVAL``` file to the directory with ```inp``` file.
1. Run ```EMCc.x``` in directory with ```inp``` and ```EIGENVAL``` files.
1. Determine real space directions from eigevectors using ```EMCcoords.pl``` script.  
*. In VASP, units are 2Pi/A.

#### 5. Acknowledgments and references
1. Mixed 2nd derivative formula: [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/).
1. Finite difference method on three-point stencil is outlined here: http://link.aip.org/link/doi/10.1063/1.2138381.

#### 6. Test cases (folders)
1. Silicon
1. Pentacene