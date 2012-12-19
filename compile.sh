#!/usr/bin/env sh

# on Cray first swap to gnu env:
#module swap PrgEnv-pgi PrgEnv-gnu
# then compile (ftn in this case is a wrapper over gfortran):
#ftn -c emc_functions.f90
#ftn emc_calc.f90 -o emc_calc.x emc_functions.o
#ftn emc_gen.f90 -o emc_gen.x emc_functions.o

# on Mac OS X Lion (and higher) or LINUX:
gfortran -llapack -lblas -c emc_functions.f90
gfortran emc_calc.f90 -llapack -lblas -o emc_calc.x emc_functions.o
gfortran emc_gen.f90 -llapack -lblas -o emc_gen.x emc_functions.o
