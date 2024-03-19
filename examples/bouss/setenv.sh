
# This is the bash file used by RJL to set the parameters needed
# to run the Bouss version of GeoClaw with MPI and OpenMP.
# Adjust as needed for your system...

export CLAW=/Users/rjl/clawpack_src/clawpack_bouss
export PYTHONPATH=$CLAW
echo CLAW is set to $CLAW

export PETSC_DIR=/Users/rjl/git/Clones/petsc
export PETSC_ARCH=arch-darwin-c-opt
export PETSC_OPTIONS="-options_file $CLAW/geoclaw/examples/bouss/petscMPIoptions"
export OMP_NUM_THREADS=6

# set prompt to show working directory and reminder bouss parameters are set:
PS1='[\W] bouss $ '
