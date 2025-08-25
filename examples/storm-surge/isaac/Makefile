
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
# If using NetCDF reading use these flags
# FFLAGS += -DNETCDF -I$(NETCDF4_DIR)/include -L$(NETCDF4_DIR)/lib
# LFLAGS += -I$(NETCDF4_DIR)/include -L$(NETCDF4_DIR)/lib -lnetcdff
FFLAGS += -DNETCDF $(shell nf-config --fflags)
LFLAGS += $(FFLAGS) $(shell nf-config --flibs) $(shell nc-config --libs)

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = 

EXCLUDE_SOURCES = 

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

RIEMANN = $(CLAW)/riemann/src

MODULES =

SOURCES = \
	$(RIEMANN)/rpn2_geoclaw.f \
	$(RIEMANN)/rpt2_geoclaw.f \
	$(RIEMANN)/geoclaw_riemann_utils.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

### DO NOT remove this line - make depends on it ###
