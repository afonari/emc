## Effective Mass Calculator Testing: fcc Si

#### 1. VASP

##### SCF
INCAR:  
```
System = fcc Si
ISTART = 0 ; ICHARG=2
ENCUT  = 300
EDIFF = 1.0E-8
ISMEAR = 0; SIGMA = 0.05;
ALGO = Fast
ISYM = 0
LCHARG = .TRUE.
LWAVE = .FALSE.
```

POSCAR:  
```
Si
1.0
        5.4306998253         0.0000000000         0.0000000000
        0.0000000000         5.4306998253         0.0000000000
        0.0000000000         0.0000000000         5.4306998253
   Si
    8
Direct
     0.000000000         0.000000000         0.000000000
     0.000000000         0.500000000         0.500000000
     0.500000000         0.500000000         0.000000000
     0.500000000         0.000000000         0.500000000
     0.750000000         0.250000000         0.750000000
     0.250000000         0.250000000         0.250000000
     0.250000000         0.750000000         0.750000000
     0.750000000         0.750000000         0.250000000
```

KPOINTS:  
```
K-Points
 0
Monkhorst Pack
 11 11 11
 0  0  0
```

OSZICAR last line: ```1 F= -.43354085E+02 E0= -.43354085E+02  d E =-.232531E-12```

##### Effective masses: Light Hole at Gamma
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
Reference value: ```0.16```.

#### 5. Acknowledgments and References
1. Silicon CIF files: [AMS](http://rruff.geo.arizona.edu/AMS/result.php?mineral=silicon)
1. Silicon effective masses: [Bart J. Van Zeghbroeck](http://ecee.colorado.edu/~bart/book/effmass.htm).
1. Effective mass of electrons in silicon: [QuantumWise tutorial](http://quantumwise.com/publications/tutorials/mini-tutorials/135-effective-mass-of-electrons-in-silicon)