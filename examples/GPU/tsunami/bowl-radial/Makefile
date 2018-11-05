
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

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

FC = pgi
USE_CUDA = TRUE
# FC = gfortran
# USE_CUDA = FALSE
 
 
# DEBUG = TRUE
DEBUG = FALSE

CC = pgi
USE_CPP = TRUE

PROFILE = TRUE
# PROFILE = FALSE 

USE_OPENMP = TRUE

VERBOSE = FALSE

# USE_DOUBLE_PRECISION = FALSE
USE_DOUBLE_PRECISION = TRUE

DIM = 2
NUM_EQUATIONS = 3
NUM_WAVES = 3
NUM_COEFFS = 1
USE_FWAVES = TRUE
USE_CAPA = FALSE

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
  $(CLAW)/riemann/src/shallow_topo_device_module.f90

SOURCES = \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f90

ifeq ($(USE_CUDA), TRUE) 
    CPP_SOURCES += stepgrid_cuda.cu
endif

# case specific dependency
$(CLAW)/riemann/src/shallow_topo_device_module.f90:$(CLAW)/amrclaw/src/2d/cudaclaw/arrayIndex.H

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

