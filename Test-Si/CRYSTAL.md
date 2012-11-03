### 1. SCF
Note:
 - P1 symmetry, don't want to mess with CRYSTAL-only trnasformation to "CONVENTIONAL" unit cell.
 - Basis set is from [Tutorial page](http://www.theochem.unito.it/crystal_tuto/mssc2008_cd/tutorials/basis_set/basis_set_tut.html).
input.d12:  
```
TEST10 - SILICON BULK - BASIS SET 6-21 MODIFIED
CRYSTAL
0 0 0
1
5.4306998253 5.4306998253 5.4306998253 90.000 90.000 90.000
8
14    0.000000000         0.000000000         0.000000000
14    0.000000000         0.500000000         0.500000000
14    0.500000000         0.500000000         0.000000000
14    0.500000000         0.000000000         0.500000000
14    0.750000000         0.250000000         0.750000000
14    0.250000000         0.250000000         0.250000000
14    0.250000000         0.750000000         0.750000000
14    0.750000000         0.750000000         0.250000000
END
14 4
2 0 6  2.  1.
2 1 6  8.  1.
2 1 2  4.  1.
0 1 1  0.  1.
0.1233392 1. 1.
99 0
END
SCFDIR
DFT
EXCHANGE
PBE
CORRELAT
PBE
END
SHRINK
8 8
TOLINTEG
8 8 8 8 16
TOLDEE
10
FMIXING
30
END
```

```== SCF ENDED - CONVERGENCE ON ENERGY      E(AU) -2.3147994980024E+03 CYCLES  12```

### 2. Effective masses: Light Hole at Gamma
inp file:
```
0.000 0.000 0.000
0.0055
1
C
```

Generating ```EIGENVAL``` file (only CRYSTAL):  
```cry-getE.pl -i ../input.out -f ../input.f9 -b 56```

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
1
C
```

Generating ```EIGENVAL``` file (only CRYSTAL):  
```cry-getE.pl -i ../input.out -f ../input.f9 -b 56```

Ouput of EMCc.x:
```
Eigensystem:
  -0.2543931563:    0.2340021460  -0.0180563778  -0.9720683941
  -0.2545028801:   -0.9569256747   0.1724609900  -0.2335603992
  -0.2545271568:   -0.1718611324  -0.9848508385  -0.0230776303
```
Reference value: ```0.25m0```.