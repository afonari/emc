### 1. SCF
Note:
 - P1 symmetry, don't want to mess with CRYSTAL-only transformation to "CONVENTIONAL" unit cell.
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
 - from output: ```== SCF ENDED - CONVERGENCE ON ENERGY      E(AU) -2.3147994980024E+03 CYCLES  12```
 - from output: ```N. OF ELECTRONS PER CELL  112```

### 2. Effective masses: Light Hole (CB-2) at Gamma
inp file:
```
0.000 0.000 0.000
0.0055
1
C
```

Generating ```EIGENVAL``` file (only CRYSTAL):  
```cry-getE.pl -i ../input.out -f ../input.f9 -b 54```

Ouput of EMCc.x:
```
 Eigensystem:
  -0.1813779485:   -0.2307390047   0.0195876899   0.9728185001
  -0.1813218118:    0.9263913328  -0.3013565081   0.2257949370
  -0.1813147635:   -0.2975879874  -0.9533103259  -0.0513888346
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