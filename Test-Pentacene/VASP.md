### 1. SCF
Check ```INCAR```, ```POSCAR``` and ```KPOINTS``` files.  
```POSCAR``` was generated from CSD: ```PENCEN```.
OSZICAR last line: ```1 F= -.49593848E+03 E0= -.49593848E+03  d E =-.541128E-13```

### 2. Effective masses: Hole at **B** ([0.375, 0.5, 0.075])
Getting reciprocal cartesian coordinates of **B** point (```T(g)*kp:``` line in output), using ```EMCcoords.pl```, ```OUTCAR``` and ```inp``` (with dummy coordinates at the first line) files are required:  
```EMCcoords.pl -k 0.375 0.5 0.075```

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
Reference value: ```0.16m0```.

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
1. Experimental and calculated effective masses: [V. Podzorov, et al.]
1. Silicon effective masses: [Bart J. Van Zeghbroeck](http://ecee.colorado.edu/~bart/book/effmass.htm).
1. Effective mass of electrons in silicon: [QuantumWise tutorial](http://quantumwise.com/publications/tutorials/mini-tutorials/135-effective-mass-of-electrons-in-silicon).
1. Silicon packing and first Brillouin zone cartoons: [S. Dhar, PhD](http://www.iue.tuwien.ac.at/phd/dhar/node18.html).