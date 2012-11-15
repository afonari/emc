#!/bin/sh

# on Mac OS X Lion (and higher) or LINUX:

gfortran -c emc_functions.f90
gfortran emc_calc.f90 -llapack -lblas -o emc_calc.x emc_functions.o
gfortran emc_gen.f90 -llapack -lblas -o emc_gen.x emc_functions.o
