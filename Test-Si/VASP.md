### 1. SCF
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
 - froum OUTCAR: ```NELECT =      32.0000```
 - OSZICAR last line: ```1 F= -.43354085E+02 E0= -.43354085E+02  d E =-.232531E-12```

### 2. Effective masses: Light Hole (CB-2) at Gamma
inp file:
```
0.000 0.000 0.000
0.01
14
V
```

Ouput of EMCc.x:
```
Eigensystem:
  -0.1615334411:    1.0000000000   0.0000000000   0.0000000000
  -0.1615334411:    0.0000000000   1.0000000000   0.0000000000
  -0.1615334411:    0.0000000000   0.0000000000   1.0000000000
```
Reference value: ```0.16m0```.

### 3. Effective masses: Heavy Hole (CB) at Gamma
```Note: unlike in the case of CRYSTAL, same EIGENVAL file can be used  to calculate effective masses for different bands, when using VASP.```

inp file:
```
0.000 0.000 0.000
0.01
16
V
```

Ouput of EMCc.x:
```
Eigensystem:
  -0.2459597491:    0.7071067812   0.5000000000  -0.5000000000
  -0.2459752609:    0.7071067812  -0.5000000000   0.5000000000
  -0.2459675047:    0.0000000000   0.7071067812   0.7071067812
```
Reference value: ```0.25m0```.