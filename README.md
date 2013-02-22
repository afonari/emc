## Effective Mass Calculator for Semiconductors

For the theory behind the concept and validation of the results please see our [arXiv paper](http://arxiv.org/abs/1302.4996).  
If you use the code (or anything else from this repo), please cite accordingly: [RIS](citation.ris) and [BibTeX](citation.bib).

### 1. Notes on theory
 1. Atomic units (a.u.) are used throughout the code: hbar = 1, energy is in Hartrees, distance is in Bohrs, mass is in the electron mass (m0).
 1. For the top of the VB (valence band) eigenvalues are negative, for the bottom of the CB (conduction band) eigenvalues are positive.
 1. In some cases, not all eigenvalues have the same sign, meaning that choosen reciprocal point is not a global minimum (maximum).
 1. Eigenvectors are directions of principal effective mass components.
 1. Effective masses (eigenvalues of the tensor) can be highly anisotropic.

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
 - **step size** in 1/Bohr. If *program identifier* is ```V``` (for *VASP*) step size will be converted to 2π/A units. At this time, in *POSCAR* scale (2nd line) should be set to **1.000**.
 - **band number**. If *CRYSTAL* is employed, band number should be set to **1**. Helper script ```cry-getE.pl``` reads-in the desired band number. For *VASP* valence band number can be obtained as ```NELECT/2``` variable from the *OUTCAR* file.
 - **program identifier** at this time can be either ```C``` (for *CRYSTAL*) or ```V``` (for *VASP*).
 - **direct lattice components** in *CRYSTAL* can be found under: ```DIRECT LATTICE VECTORS COMPON. (A.U.)```. In *VASP* under: ```direct lattice vectors```. Program will deal with units itself.
 - Please remove comments from ```inp``` file, as they are not currently supported.

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