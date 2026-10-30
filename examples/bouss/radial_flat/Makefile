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


EXE = $(PWD)/xgeoclaw
SETRUN_FILE = setrun.py           # File containing function to make data
OUTDIR = _output                  # Directory for output
SETPLOT_FILE = setplot.py         # File containing function to set plots
PLOTDIR = _plots                  # Directory for plots


FFLAGS ?= -O  -gno-strict-dwarf -fbounds-check  -fopenmp \
      -std=legacy -ffpe-trap='invalid,overflow,zero'


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

