
# Makefile for the topo descriptor-crop regression case.
# Builds a standard single-layer GeoClaw executable; identical to the
# met_forcing Makefile, since what is under test is the topo reader, not the
# forcing.  NetCDF is required (the descriptor path is topo_type=4 only).
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

FFLAGS ?=

# NetCDF support is opt-in so the parametric (Holland) build stays dependency
# free.  The gridded netCDF aux case builds with:
#   make new USE_NETCDF=1 NETCDF_FFLAGS="$(nf-config --fflags)" \
#            NETCDF_LFLAGS="$(nf-config --flibs)"
# (the test injects these automatically via build_executable make_vars).
ifeq ($(USE_NETCDF),1)
FFLAGS += -DNETCDF $(NETCDF_FFLAGS)
LFLAGS += $(FFLAGS) $(NETCDF_LFLAGS)
endif

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

RIEMANN = $(CLAW)/riemann/src

MODULES = \

SOURCES = \
  $(RIEMANN)/rpn2_geoclaw.f \
  $(RIEMANN)/rpt2_geoclaw.f \
  $(RIEMANN)/geoclaw_riemann_utils.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

### DO NOT remove this line - make depends on it ###
