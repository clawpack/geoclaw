# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use

# These environment variables need to be set properly, usually in your shell:
#PETSC_DIR=/Users/rjl/git/Clones/petsc
#PETSC_ARCH=arch-darwin-c-opt
#PETSC_OPTIONS=-options_file $(PETSC_DIR)/petscMPIoptions
PETSC_DIR=/Users/berger/petsc
PETSC_ARCH=arch-mpi-opt

EXE = $(PWD)/xgeoclaw
#RUNEXE="${EXE}"
RUNEXE="${PETSC_ARCH}/bin/mpiexec -n 6 ${EXE}"
SETRUN_FILE = setrun.py           # File containing function to make data
OUTDIR = _output                  # Directory for output
SETPLOT_FILE = setplot.py         # File containing function to set plots
PLOTDIR = _plots          # Directory for plots


# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FC = gfortran
PETSC_INCLUDE = $(PETSC_DIR)/include $(PETSC_DIR)/$(PETSC_ARCH)/include         
INCLUDE += $(PETSC_INCLUDE)
# line below returns extra -I -I...
#PETSC_INCLUDE = $(shell PKG_CONFIG_PATH=$(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig pkg-config --cflags PETSc)
FFLAGS ?= -O  -gno-strict-dwarf -fbounds-check  -fopenmp \
      -std=legacy -ffpe-trap='invalid,overflow,zero'
FFLAGS += -DHAVE_PETSC -ffree-line-length-none
#LFLAGS += -Wl,-rpath,/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -Wl,-rpath,/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -L/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -Wl,-rpath,/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -L/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -Wl,-rpath,/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12 -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12 -Wl,-rpath,/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc -Wl,-rpath,/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current -lpetsc -llapack -lblas -lc++ -ldl -lmpifort -lmpi -lpmpi -lgfortran -lemutls_w -lquadmath -lc++ -ldl
#LFLAGS +=   -Wl,-rpath,$(PETSC_ARCH)/lib -L$(PETSC_ARCH)/lib -Wl,-rpath,/opt/X11/lib -L/opt/X11/lib -Wl,-rpath,/opt/intel/compilers_and_libraries_2020.2.258/mac/ipp/lib -L/opt/intel/compilers_and_libraries_2020.2.258/mac/ipp/lib -Wl,-rpath,/opt/intel/compilers_and_libraries_2020.2.258/mac/compiler/lib -L/opt/intel/compilers_and_libraries_2020.2.258/mac/compiler/lib -Wl,-rpath,/opt/intel/compilers_and_libraries_2020.2.258/mac/mkl/lib -L/opt/intel/compilers_and_libraries_2020.2.258/mac/mkl/lib -Wl,-rpath,/opt/intel/compilers_and_libraries_2020.2.258/mac/tbb/lib -L/opt/intel/compilers_and_libraries_2020.2.258/mac/tbb/lib -Wl,-rpath,/opt/intel/compilers_and_libraries_2020.2.258/mac/daal/lib -L/opt/intel/compilers_and_libraries_2020.2.258/mac/daal/lib -Wl,-rpath,/usr/local/Cellar/gcc@8/8.4.0_2/lib/gcc/8/gcc/x86_64-apple-darwin19/8.4.0 -lpetsc -llapack -lblas -lX11 -ldl -lgfortran -lquadmath -lm -lquadmath -ldl
PETSC_LFLAGS = $(shell PKG_CONFIG_PATH=$(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig pkg-config --libs-only-L --libs-only-l PETSc)
LFLAGS += $(PETSC_LFLAGS) -fopenmp

# ---------------------------------
# package sources for this program:
# ---------------------------------

# BOUSSLIB contains library code for Boussinesq solvers,
# AMRLIB and GEOLIB are standard libraries, defined in case you need to
# exclude some modele or source file from one of these.
BOUSSLIB = $(CLAW)/geoclaw/src/2d/bouss
AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow

# See this Makefile for the list of library routines used:
include $(BOUSSLIB)/Makefile.bouss

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------


MODULES = \

SOURCES = \
  ./qinit.f90 \



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

