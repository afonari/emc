## Effective Mass Calculator

#### 1. Theory

```
Note: Will work in atomic units (a.u.), where: hbar = 1, energy is in Hartrees,
distance is in Bohrs, mass is in electron mass (m0).
```

Effective mass tensor is defined as:  
![Effective Mass Tensor](https://raw.github.com/alexandr-fonari/emc/master/p_ms.gif)

Energy gradient tensor (symmetric) is defined in 3D as follows:  
![Energy Gradient Tensor](https://raw.github.com/alexandr-fonari/emc/master/p_et.png)  
where ```x*, y*, z*``` are cartesian reciprocal directions (see ```EMCcoords.py``` for reference).

Second order derivatives are estimated with ```O(h^4)``` error:  
![2nd Derivative](https://raw.github.com/alexandr-fonari/emc/master/p_2ndd.png)

Mixed second derivatives are estimated also with ```O(h^4)``` error:  
![Mixed 2nd Derivative](http://www.holoborodko.com/pavel/wp-content/ql-cache/quicklatex.com-ead43440eddb0f8db2cc36a1df79c547_l3.svg)

After matrix is built, it is inversed (```DGETRF+DGETRI``` [[SO]](http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c)) and diagonalized (```DGEEV```).

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
1. Make directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for **Y** point with ```dk=0.01```.
1. Copy ```KPOINTS``` file generated with ```EMCg.x``` into it.
1. Run from freshly created directory (```emH-00-50-00-d01```) as ```cry-getE.pl -i ../input.out -f ../input.f9 -b 131```.
1. Copy ```EIGENVAL``` generated file to the directory with ```inp``` file.
1. Run ```EMCc.x``` in directory with ```inp``` and ```EIGENVAL``` files.

#### 4. How to run with VASP
1. Run SCF (e.g. ```ICHARG=2``` in INCAR).
1. Make directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```.
1. Copy ```KPOINTS``` file generated with ```EMCg.x``` into it.
1. Run non-SCF calculation ( ```ICHARG=11``` in INCAR) from freshly created directory (```emH-00-50-00-d01```)
1. Copy ```EIGENVAL``` generated file to the directory with ```inp``` file.
1. Run ```EMCc.x``` in directory with ```inp``` and ```EIGENVAL``` files.

#### 5. Acknowledgments and References
1. Mixed 2nd derivative formula: [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/)