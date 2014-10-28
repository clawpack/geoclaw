
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = classic                  # Clawpack package to use
EXE = xclaw                         # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=

# ---------------------------------
# List of sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/1d/shallow

MODULES = \
  $(GEOLIB)/gauges_module.f90 \
  $(GEOLIB)/geoclaw_module.f90

SOURCES = \
  qinit.f \
  setaux.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \
  $(GEOLIB)/rp1_geoclaw.f90 \
  $(GEOLIB)/setprob.f90 \
  $(GEOLIB)/b4step1.f90 \
  $(GEOLIB)/bc1.f \
  $(GEOLIB)/out1.f \
  $(GEOLIB)/src1.f90 \
  $(CLAW)/classic/src/1d/driver.f90 \
  $(CLAW)/classic/src/1d/inlinelimiter.f90 \
  $(CLAW)/classic/src/1d/step1.f90 \
  $(CLAW)/classic/src/1d/claw1ez.f \
  $(CLAW)/classic/src/1d/claw1.f \
  $(CLAW)/classic/src/1d/copyq1.f \
  $(CLAW)/classic/src/1d/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

