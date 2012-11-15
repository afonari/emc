#!/bin/sh

# on Mac OS X Lion (and higher) or LINUX:

gfortran -c emc_functions.f90
# gfortran EMCc.f90 -llapack -lblas -o EMCc.x emc_functions.o
gfortran EMCg.f90 -llapack -lblas -o EMCg.x emc_functions.o
