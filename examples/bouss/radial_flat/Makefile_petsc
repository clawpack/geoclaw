# Makefile for Clawpack code in this directory.

# Version for Boussinesq solvers in GeoClaw, which requires PETSc and MPI
# See https://www.clawpack.org/bouss2d.html
#
# Execute this command to check that enviornment variables set properly:
#    make check
# before compiling or executing the code!
# Then do `make new` if using the Boussinesq version for the first time.

# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use

# Several environment variables are expected in this Makefile,
# as set for example in the bash script
#     $CLAW/geoclaw/examples/bouss/setenv.sh
# Some of these could instead be set here explicitly.

# But PETSC_OPTIONS must be set as an environment variable!
# The file petscMPIoptions sets parameters for MPI and the iterative solver:
# e.g. in bash:
#   export PETSC_OPTIONS="-options_file $CLAW/geoclaw/examples/bouss/petscMPIoptions"

ifndef PETSC_DIR
  $(error PETSC_DIR not set)
endif

ifndef PETSC_OPTIONS
  PETSC_OPTIONS=MISSING
  $(error PETSC_OPTIONS must be declared as environment variable)
endif

ifndef CLAW_MPIEXEC
  CLAW_MPIEXEC=MISSING
  $(error CLAW_MPIEXEC must be declared as environment variable)
endif

ifndef CLAW_MPIFC
  CLAW_MPIFC=MISSING
  $(error Fortran compiler CLAW_MPIFC must be declared as environment variable)
endif

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC = ${CLAW_MPIFC}

# How many MPI processes to use:
BOUSS_MPI_PROCS ?= 6

EXE = $(PWD)/xgeoclaw
RUNEXE="${CLAW_MPIEXEC} -n ${BOUSS_MPI_PROCS}"
SETRUN_FILE = setrun.py           # File containing function to make data
OUTDIR = _output                  # Directory for output
SETPLOT_FILE = setplot.py         # File containing function to set plots
PLOTDIR = _plots                  # Directory for plots


# Some compiler flags below are needed for PETSc
PETSC_INCLUDE = $(PETSC_DIR)/include $(PETSC_DIR)/$(PETSC_ARCH)/include
INCLUDE += $(PETSC_INCLUDE)
PETSC_LFLAGS = $(shell PKG_CONFIG_PATH=$(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig pkg-config --libs-only-L --libs-only-l PETSc)

FFLAGS ?= -O  -gno-strict-dwarf -fbounds-check  -fopenmp \
      -std=legacy -ffpe-trap='invalid,overflow,zero'
FFLAGS += -DHAVE_PETSC -ffree-line-length-none
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

.PHONY: check

check:
	@echo ===================
	@echo CLAW = $(CLAW)
	@echo OMP_NUM_THREADS = $(OMP_NUM_THREADS)
	@echo BOUSS_MPI_PROCS = $(BOUSS_MPI_PROCS)
	@env | grep PETSC_OPTIONS
	@echo PETSC_DIR = $(PETSC_DIR)
	@echo PETSC_ARCH = $(PETSC_ARCH)
	@echo INCLUDE = $(INCLUDE)
	@echo CLAW_MPIEXEC = $(CLAW_MPIEXEC)
	@echo RUNEXE = $(RUNEXE)
	@echo EXE = $(EXE)
	@echo CLAW_MPIFC = $(CLAW_MPIFC)
	@echo FC = $(FC)
	@echo FFLAGS = $(FFLAGS)
	@echo LFLAGS = $(LFLAGS)
	@echo OUTDIR = $(OUTDIR)
	@echo PLOTDIR = $(PLOTDIR)
	@echo ===================
