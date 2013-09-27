---
layout: default
title: Installing NWChem 6.3 on Cray XE
---

## {{ page.title }}

### Compiling NWChem
 1. Download NWChem ([version](http://www.nwchem-sw.org/images/Nwchem-6.3.revision1-src.2013-05-28.tar.gz) used in this post) in the home directory
 1. Unpack it:  
 `tar -zxvf Nwchem-6.3.revision1-src.2013-05-28.tar.gz`
 1. Change directory to the NWChem directory:  
 `cd nwchem-6.3-src.2013-05-28/`
 1. Define environment variables (see below contents of the vars.sh file): `source vars.sh`
 1. Check that correct modules are loaded (see below)
 1. Change directory to *src*: `cd src/`
 1. Run: `make 64_to_32`
 1. Run (*NWCHEM_CONFIG* file will be created automatically): `gmake >& make.log`

### Running NWChem

Example of PBS script:

```bash
#!/bin/bash
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=00:10:00
#PBS -q debug
#PBS -j oe
#PBS -N jobtitle
#PBS -m abe
#PBS -M firstame.lastname@email.com

# nwchem vars
export NWCHEM_TOP=/home/user/nwchem-6.3-src.2013-05-28
export NWCHEM_TARGET=LINUX64
export NWCHEM_BASIS_LIBRARY=$NWCHEM_TOP/src/basis/libraries/
export SCRATCH_DIR=/scratch/user
export MPI_DSM_DISTRIBUTE=1

# hugepages vars
export HUGETLB_MORECORE=yes
export HUGETLB_ELFMAP=W
export HUGETLB_FORCE_ELFMAP=yes+

cd $PBS_O_WORKDIR
nwchem_exe=${NWCHEM_TOP}/bin/${NWCHEM_TARGET}/nwchem

job_name=ethylene-dimer

aprun -B $nwchem_exe $job_name.nw > $job_name.out
```

### Some helpful reference data

#### Cray machine

Compilation variables:

```bash
export NWCHEM_TOP=`pwd`
export NWCHEM_MODULES="all"
export NWCHEM_TARGET=LINUX64
export TARGET=LINUX64
export USE_NOFSCHECK=TRUE
export LARGE_FILES=y

export USE_MPI=y
export USE_MPIF=y
export USE_MPIF4=y
export USE_SCALAPACK=y
export USE_64TO32=y
export LIBMPI=" "
export ARMCI_NETWORK=DMAPP
export ONESIDED_USE_UDREG=1

export USE_SCALAPACK=y
export SCALAPACK=-lsci_pgi
export SCALAPACK_LIB=-lsci_pgi

export HAS_BLAS=yes
export BLASOPT=-lsci_pgi
export BLAS_OPT=-lsci_pgi
export BLAS_LIB=-lsci_pgi

export LAPACK_LIB=-lsci_pgi

export FC=ftn
export CC=cc
```

Loaded modules:

```
$ module list
[...]
PrgEnv-pgi/4.0.46
onesided/1.5.0
craype-hugepages64M
[...]
```

### Notes

 * Instead of *all*, other strings can be used in `NWCHEM_MODULES` keyword to compile only desired modules (*e.g.* `qm` for quantum). See `INSTALL` file for a complete description.
 * See this [thread](https://groups.google.com/forum/#!topic/hpctools/ZvHgljFkYWg) for `ARMCI_NETWORK` type discussion for Cray XE machines.

#### Generic Linux cluster

Compilation variables:

```bash
export NWCHEM_TOP=`pwd`
export NWCHEM_MODULES="all"
export USE_NOFSCHECK=TRUE
export NWCHEM_TARGET=LINUX64

export MPI_HOME=/usr/local/packages/openmpi/1.5.4/intel-11.1.059/
export MKLROOT=/usr/local/packages/intel/mkl/10.3/lib/intel64

export USE_MPI=y
export USE_MPIF=y
export USE_MPIF4=y

export MPI_INCLUDE=${MPI_HOME}/include
export MPI_LIB=${MPI_HOME}/lib
export LIBMPI="-lmpi_f90 -lmpi_f77 -lmpi -ldl -Wl,--export-dynamic -lnsl -lutil"

export BLASOPT=" -L${MKLROOT} -Wl,--start-group -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group"
export BLAS_LIB=${BLASOPT}
export BLAS_SIZE=8

export ARMCI_NETWORK=MPI-TS
export FC=ifort
```