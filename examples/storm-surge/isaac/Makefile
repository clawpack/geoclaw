
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
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots


# Environment variable FC should be set to fortran compiler, e.g. gfortran
FFLAGS += -Wall -pedantic -fcheck=all -fbacktrace -fno-f2c -O0 -g -DNETCDF -lnetcdff -I/Users/catherinej/homebrew/Cellar/netcdf/4.8.1/include -L/Users/catherinej/homebrew/Cellar/netcdf/4.8.1/
 # Need to add LFLAGS here explicitly otherwise will not be linked
LFLAGS += -lnetcdff -L/Users/catherinej/homebrew/Cellar/netcdf/4.8.1/lib
# ---------------------------------
# package sources for this program:
# ---------------------------------

#AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

NEWAMR = ../../new_fortran_amrclaw
NEWGEO = ../../new_fortran_geoclaw
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

RIEMANN = $(CLAW)/riemann/src

MODULES = \


SOURCES = \
  $(RIEMANN)/rpn2_geoclaw.f \
  $(RIEMANN)/rpt2_geoclaw.f \
  $(RIEMANN)/geoclaw_riemann_utils.f \


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
# Construct the topography data

# Construct the topography data
.PHONY: all

all:
	$(MAKE) .plots
	$(MAKE) .htmls## DO NOT remove this line - make depends on it ###
