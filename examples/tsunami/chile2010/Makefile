# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Define GEOCLAW_VEC as "1" to use the code with vectorized Riemann solvers.
# The following compiler flags are suggested for vectorization (Intel Compiler):
#   -O2 -fopenmp -align array64byte -xHost


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= 

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

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

ifneq ($(GEOCLAW_VEC),1)
  SOURCES = \
    $(CLAW)/riemann/src/rpn2_geoclaw.f \
    $(CLAW)/riemann/src/rpt2_geoclaw.f \
    $(CLAW)/riemann/src/geoclaw_riemann_utils.f 
else
  EXCLUDE_SOURCES = \
    $(GEOLIB)/qad.f \
    $(GEOLIB)/flux2fw.f \
    $(GEOLIB)/step2.f \
    $(AMRLIB)/inlinelimiter.f

  SOURCES = \
    $(CLAW)/riemann/src/rpn2_geoclaw_vec.f90 \
    $(CLAW)/riemann/src/rpt2_geoclaw_vec.f90 \
    $(CLAW)/riemann/src/geoclaw_riemann_utils_vec.f90 \
    $(GEOLIB)/qad_vec.f \
    $(GEOLIB)/flux2fw_vec.f \
    $(GEOLIB)/step2_vec.f90 \
    $(GEOLIB)/inlinelimiter_vec.f
endif


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

