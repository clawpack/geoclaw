
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
EXE = xgeo                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS = -llapack -lblas -O2


# ---------------------------------
# package sources for this program:
# ---------------------------------

BOUSSLIB = $(CLAW)/geoclaw/src/1d_classic/bouss
include $(BOUSSLIB)/Makefile.geoclaw_1d_classic_bouss

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
  bc1.f \
  setprob.f90 \
  #$(BOUSSLIB)/lapack_tridiag.f \  # alternative to linking -lapack -blas


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all

topo:
	python make_celledges.py

all: 
	$(MAKE) topo
	$(MAKE) .plots
	$(MAKE) .htmls

