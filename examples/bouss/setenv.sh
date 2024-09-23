
# This is a sample bash script to set the parameters needed
# to run the Bouss version of GeoClaw with MPI and OpenMP.
# Adjust as needed for your system...

# You also need to set CLAW, FC, and perhaps PYTHONPATH

# For more information, see
#   https://www.clawpack.org/bouss2d.html
#   https://www.clawpack.org/setenv.html

export PETSC_DIR=/full/path/to/petsc
export PETSC_ARCH=arch-darwin-c-opt
export PETSC_OPTIONS="-options_file $CLAW/geoclaw/examples/bouss/petscMPIoptions"
export OMP_NUM_THREADS=6
export BOUSS_MPI_PROCS=6  # only used in Clawpack Boussinesq example

