### 1. SCF
Check ```INCAR```, ```POSCAR``` and ```KPOINTS``` files. ```POSCAR``` was generated from CSD: ```PENCEN```.  
OSZICAR last line: ```1 F= -.49593848E+03 E0= -.49593848E+03  d E =-.541128E-13```

### 2. Effective masses: Hole at **B** ([0.375, 0.5, 0.075])
Converted to reciprocal cartesian coordinates **B** point (look for ```T(g)*kp:``` line in the output of ```EMCcoords.pl```), using ```EMCcoords.pl```, ```OUTCAR``` and ```inp``` (with dummy coordinates at the first line) files are required:
```
EMCcoords.pl -k 0.375 0.5 0.075
```

inp file:
```
0.0474684  0.0792446  0.0409311
0.01
102
V
```

Ouput of EMCc.x:
```
 Eigensystem:
 -27.4579649685:   -0.1438683390  -0.7905320511  -0.5952822668
  -2.8049633222:   -0.9880379472   0.1484981479   0.0415850334
  -5.0732884101:   -0.0555240123  -0.5941442386   0.8024398468
```

Converted eigenvectors of the matrix to the real space directions (using ```EMCcoords.pl```, look for ```Norm(g*kp)``` line):
```
m0 =  -2.80: -0.979a + 0.204b + 0.023c
m0 =  -5.07:  0.404a - 0.726b + 0.557c
m0 = -27.46: -0.259a - 0.932b - 0.255c
```

Reference (calculated) values from [1]:
```

```

### 3. Effective masses: Heavy Hole at Gamma
inp file:
```
0.000 0.000 0.000
0.0055
16
V
```

Ouput of EMCc.x:
```
Eigensystem:
  -0.2440325848:    0.7071067812  -0.7071067812   0.0000000000
  -0.2440345045:    0.7071067812   0.7071067812   0.0000000000
  -0.2440335447:    0.0000000000   0.0000000000   1.0000000000
```
Reference value: ```0.25m0```.

### Acknowledgments and References
1. K. Doi, *et al.*, *J. Appl. Phys*, **98**, 113709 (2005): http://dx.doi.org/10.1063/1.2138381
1. Silicon effective masses: [Bart J. Van Zeghbroeck](http://ecee.colorado.edu/~bart/book/effmass.htm).
1. Effective mass of electrons in silicon: [QuantumWise tutorial](http://quantumwise.com/publications/tutorials/mini-tutorials/135-effective-mass-of-electrons-in-silicon).
1. Silicon packing and first Brillouin zone cartoons: [S. Dhar, PhD](http://www.iue.tuwien.ac.at/phd/dhar/node18.html).