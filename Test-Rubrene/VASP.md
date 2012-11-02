### 1. SCF
Check ```INCAR```, ```POSCAR``` and ```KPOINTS``` files.  
```POSCAR``` was generated from CSD: ```QQQCIG01```.
OSZICAR last line: ```1 F= -.19021982E+04 E0= -.19021982E+04  d E =-.254312E-11```

### 2. Effective masses: Light Hole at Gamma
inp file:
```
0.000 0.000 0.000
0.0055
14
V
```

Ouput of EMCc.x:
```
Eigensystem:
  -0.1620290504:    1.0000000000   0.0000000000   0.0000000000
  -0.1620290504:    0.0000000000   1.0000000000   0.0000000000
  -0.1620290504:    0.0000000000   0.0000000000   1.0000000000
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