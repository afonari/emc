## Effective Mass Calculator
### 1. Theory

```
Note: Atomic units (a.u.) will be used, hbar = 1, energy is in Hartrees,
distance is in Bohrs, mass is in electron mass (m0).
```

The average velocity of an electron for a certain value of *n* and *k* is:  
![Mean velocity](https://raw.github.com/alexandr-fonari/emc/master/pics/p_v.png)  
thus acceleration can be written as:  
![Acceleration](https://raw.github.com/alexandr-fonari/emc/master/pics/p_a.png)

In the semiclassical picture, equation of motion of an electron in the presence of the electric ( **E** ) and magnetic ( **H** ) fields is:  
![Equation of motion](https://raw.github.com/alexandr-fonari/emc/master/pics/p_e_m.png)  
in the case of: ```H = 0```:  
![Equation of motion](https://raw.github.com/alexandr-fonari/emc/master/pics/p_f.png)  
getting back to acceleration:  
![Acceleration](https://raw.github.com/alexandr-fonari/emc/master/pics/p_a2.png)  
from here, inverse effective mass tensor (9 comoponents) is expressed through the tensor of second derivatives of energy with respect to reciprocal wavevector:  
![Inverse EM tensor](https://raw.github.com/alexandr-fonari/emc/master/pics/p_1o_m.png)  
Note that this tensor is symmetric (can be diagonalized with [DSYEV](http://netlib.org/lapack/double/dsyev.f)):  
![Tensor](https://raw.github.com/alexandr-fonari/emc/master/pics/p_tensor.png)  
where ```x*, y*, z*``` are reciprocal directions.

At a saddle point (e.g. band maximum/minimum) components of the effective mass are inverse of eigenvalues of the tensor:  
![Eigenvalues](https://raw.github.com/alexandr-fonari/emc/master/pics/p_ev.png)  
where ```αi``` are eigevalues of the ```d²E/dk²``` tensor.

**Notes**:
 1. For the top of the band, eigenvalues are negative, for the bottom of the band, eigenvalues are positive.
 1. In some cases, not all eigenvalues have the same sign, meaning that choosen reciprocal point is not a global minimum (maximum).
 1. Eigenvectors are directions of principal effective mass components.
 1. Note that components of the effective mass tensor can be highly anisotropic.

#### 1.1 Numerical Differentiation
Derivatives are estimated using finite-difference method on a five-point stencil.  
Second order derivatives are estimated with ```O(h^4)``` error:  
![2nd Derivative](https://raw.github.com/alexandr-fonari/emc/master/pics/p_2ndd.png)  
Mixed second derivatives are estimated also with ```O(h^4)``` error:  
![Mixed 2nd Derivative](http://www.holoborodko.com/pavel/wp-content/ql-cache/quicklatex.com-ead43440eddb0f8db2cc36a1df79c547_l3.svg)

### 2. Required input files
```inp``` has the following form:  
```
0.000 0.000 0.000                       ! K-POINT in reciprocal space (3 floats)
0.01                                    ! step size (1 float)
81                                      ! band number, (1 integer)
V                                       ! program indentifier (1 char)
6.291999817  0.000000000  0.000000000   ! direct lattice vectors
0.755765092  7.652872670  0.000000000   ! direct lattice vectors
0.462692761  3.245907103 14.032346772   ! direct lattice vectors
```
 - **K-POINT** coordinates in reciprocal space of a band maximum for holes and band minimum for electrons.
 - **step size** in 1/Bohr. If *program identifier* is ```V``` (for *VASP*) step size will be converted to 2π/A units. At this time, in *POSCAR* scale (2nd line) should be set to **1.**.
 - **band number**. If *CRYSTAL* is employed, band number should be set to **1**. Helper script ```cry-getE.pl``` reads in desired band number. For *VASP* valence band number can be obtained as ```NELECT/2``` variable from the *OUTCAR* file.
 - **program identifier** at this time can be either ```C``` (for *CRYSTAL*) or ```V``` (for *VASP*).
 - **direct lattice components** in *CRYSTAL* can be found under: ```DIRECT LATTICE VECTORS COMPON. (A.U.)```. In *VASP* under: ```direct lattice vectors```. Program will deal with units.
 - Please remove comments from ```inp``` file, as they are not currently supported.

### 3. Usage
#### 3.1 CRYSTAL
1. Run SCF.
1. Create a directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```.
1. Make ```inp``` file with the desired characteristics in the newly created directory.
1. Run ```emc_gen.x``` to obtain ```KPOINTS``` file.
1. Run ```cry-getE.pl -f ../input.f9 -b 131``` to obtain ```EIGENVAL``` file.
 - **cry-getE.pl** uses the following parameters:
   * SCF f9 file in ```-f``` flag;
   * band number in ```-b``` flag.
   * Note that ```runprop09``` needs to be in the current ```$PATH```, otherwise script will quit.
1. Run ```emc_calc.x``` to obtain effective masses and directions. Look for ```emc_calc.log``` file.

#### 3.2 VASP
1. Run SCF (e.g. ```ICHARG=2``` in INCAR).
1. Create a directory, e.g. ```emH-00-50-00-d01```, meaning we are calculating effective mass for holes (VB) **Y** point with ```dk=0.01```.
1. Make ```inp``` file with the desired characteristics in the newly created directory.
1. Run ```emc_gen.x``` to obtain ```KPOINTS``` file.
1. Run a non-SCF calculation (```ICHARG=11``` in INCAR). Don't forget to copy ```CHGCAR``` file from SCF folder.
1. Run ```emc_calc.x``` to obtain effective masses and directions. Look for ```emc_calc.log``` file.

### 4. Acknowledgments and references
1. Mixed 2nd derivative formula: [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/).
1. Finite difference method on three-point stencil is outlined here: http://link.aip.org/link/doi/10.1063/1.2138381.

### 5. Test cases (folders)
1. Silicon
1. Pentacene