
# This is a sample bash script to set the parameters needed
# to run the Bouss version of GeoClaw with MPI and OpenMP.
# Adjust as needed for your system...

# For more information on installing PETSc and setting these environment
# variables properly, see $CLAW/geoclaw/examples/bouss/README.txt and
#   https://www.clawpack.org/bouss2d.html

# You also need to set CLAW, and perhaps PYTHONPATH, see:
#   https://www.clawpack.org/setenv.html

# export CLAW=full/path/to/clawpack
echo CLAW is set to $CLAW

# path to PETSc installation:
export PETSC_DIR=/full/path/to/petsc

# PETSC_ARCH is only needed if PETSc is installed inside the PETSc directory.
export PETSC_ARCH=arch-darwin-c-opt
# For PETSc installs by conda or package managers, it should not be set.
#export PETSC_ARCH=

# PETSC_OPTIONS should point to the options file needed to specify
# many parameters controlling which solver, preconditioner, etc are used:
export PETSC_OPTIONS="-options_file $CLAW/geoclaw/examples/bouss/petscMPIoptions"

# number of OpenMP threads to use for explicit AMR (as usual in GeoClaw):
export OMP_NUM_THREADS=6

# number of MPI processes to use for solving linear systems with PETSc:
export BOUSS_MPI_PROCS=6

# the mpiexec command to use to run the executable:
export CLAW_MPIEXEC=mpiexec
# set CLAW_MPIEXEC to mpiexec only if this command is defined in your shell,
# e.g. to use some version of MPI was installed outside of PETSc.
# Or set to the full path to this command, e.g. for the PETSc version:
#export CLAW_MPIEXEC=$PETSC_DIR/$PETSC_ARCH/bin/mpiexec  # requires PETSC_ARCH

# the proper Fortran compiler to use for MPI code:
# e.g. mpif90 if that is defined in your shell, or gfortran *might* work.
# This will over-rule any FC environment variable.
export CLAW_MPIFC=mpif90
