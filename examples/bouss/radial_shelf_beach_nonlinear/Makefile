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

#PETSC_DIR=/Users/rjl/git/Clones/petsc
#PETSC_ARCH=$(PETSC_DIR)/arch-darwin-c-opt
#PETSC_OPTIONS=-options_file $(PETSC_DIR)/petscMPIoptions # from MJB

EXE = $(PWD)/xgeoclaw
RUNEXE="${PETSC_ARCH}/bin/mpiexec -n 6 ${EXE}"
SETRUN_FILE = setrun.py           # File containing function to make data
OUTDIR = _output                  # Directory for output
SETPLOT_FILE = setplot.py         # File containing function to set plots
PLOTDIR = _plots                  # Directory for plots


# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FC = gfortran
INCLUDE += $(PETSC_DIR)/include
INCLUDE += $(PETSC_ARCH)/include         
FFLAGS =  -DHAVE_PETSC  -gno-strict-dwarf -fbounds-check -fopenmp \
      -ffree-line-length-none -std=legacy
LFLAGS ?=
LFLAGS += -fopenmp
LFLAGS += -Wl,-rpath,/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -Wl,-rpath,/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -L/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -Wl,-rpath,/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -L/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -Wl,-rpath,/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12 -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12 -Wl,-rpath,/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc -Wl,-rpath,/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current -lpetsc -llapack -lblas -lc++ -ldl -lmpifort -lmpi -lpmpi -lgfortran -lemutls_w -lquadmath -lc++ -ldl



# ---------------------------------
# package sources for this program:
# ---------------------------------

#GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
#BOUSSLIB = $(CLAW)/geoclaw/src/2d/bouss

# for testing:
GEOLIB = /Users/rjl/clawpack_src/clawpack_bouss/geoclaw/src/2d/shallow
BOUSSLIB = /Users/rjl/clawpack_src/clawpack_bouss/geoclaw/src/2d/bouss

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
  ./flag2refine2.f90 \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all
topo:
	python maketopo.py

all: 
	$(MAKE) topo
	$(MAKE) .plots
	$(MAKE) .htmls

run:
	$(RUNEXE)
