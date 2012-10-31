#!/bin/sh

# on Mac OS X Lion (and higher) or LINUX:

gfortran EMCg.f90 -o EMCg.x
gfortran EMCc.f90 -llapack -lblas -o EMCc.x