LAPACKLIBS=-llapack -lblas #  [for Mac os] comment this line if you are not using LAPACK.
                             # On systems with a FORTRAN (not f2c'ed) version of LAPACK, -lf2c is
                             # not necessary; on others, -lf2c is equivalent to -lF77 -lI77

#LAPACKLIBS=-L/usr/local/atlas/lib -llapack -lcblas -lf77blas -latlas -lf2c # This works with   the ATLAS updated lapack and Linux_P4SSE2
                                    # from   http://www.netlib.org/atlas/archives/linux/

#LAPACKLIBS=-llapack -lgoto2 -lpthread -lf2c # This works with GotoBLAS
                                         # from http://www.tacc.utexas.edu/research-development /tacc-projects/

#LAPACKLIBS=-L/opt/intel/mkl/8.0.1/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lf2c # This works with MKL 8.0.1 from
                   # http://www.intel.com/cd/software/products/asmo-na/eng/perflib/mkl/index.htm

#LAPACKLIBS=-L/opt/intel/Compiler/11.1/046/mkl/lib/em64t -lmkl_lapack -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread

LIBS=$(LAPACKLIBS)

#all of these Fortran compilers work fine for cuda/cula programs
#note the option for the g95 compiler
FOR=gfortran

all: emc_functions.o emc_gen emc_calc

emc_functions.o: emc_functions.f90
	$(FOR) $(LIBS) -c emc_functions.f90

emc_gen: emc_gen.f90 emc_functions.o
	$(FOR) emc_gen.f90 $(LIBS) -o emc_gen emc_functions.o

emc_calc: emc_calc.f90 emc_functions.o
	$(FOR) emc_calc.f90 $(LIBS) -o emc_calc emc_functions.o

clean:
	@rm -f *.o *.mod emc_calc emc_gen

#end