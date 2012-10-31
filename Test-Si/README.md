## Effective Mass Calculator Testing: fcc Si

#### 1. VASP

##### SCF
INCAR:  
```
System = fcc Si
ISTART = 0 ; ICHARG=2
ENCUT  =    240
ISMEAR = 0; SIGMA = 0.1;
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

#### 5. Acknowledgments and References
1. Mixed 2nd derivative formula: [Pavel Holoborodko](http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/)