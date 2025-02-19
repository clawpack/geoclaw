
# This is a sample bash script to set the parameters needed
# to run the Bouss version of GeoClaw with MPI and OpenMP.
# Adjust as needed for your system...

# For more information, see
#   https://www.clawpack.org/bouss2d.html
#   https://www.clawpack.org/setenv.html

# You also need to set CLAW and perhaps PYTHONPATH, see:
#   https://www.clawpack.org/setenv.html

echo CLAW is set to $CLAW

# path to PETSc installation:
export PETSC_DIR=/full/path/to/petsc  # NEED TO FIX!

# PETSC_ARCH is only needed if PETSc is installed inside the PETSc directory.
# For PETSc installs by conda or package managers, it should not be set.
#export PETSC_ARCH=
export PETSC_ARCH=arch-darwin-c-opt  # NEED TO FIX!

# You may want to use a different version of petscMPIoptions
# This setting uses the version in this directory:
export PETSC_OPTIONS="-options_file $CLAW/geoclaw/examples/bouss/petscMPIoptions"

export OMP_NUM_THREADS=6
export BOUSS_MPI_PROCS=6

# CLAW_MPIEXEC should be set to the command used to execute MPI code:
export CLAW_MPIEXEC=mpiexec
# set CLAW_MPIEXEC to mpiexec only if this command is defined in your shell,
# e.g. to use some version of MPI was installed outside of PETSc.
# Or set to the full path to this command, e.g. for the PETSc version:
#export CLAW_MPIEXEC=$PETSC_DIR/$PETSC_ARCH/bin/mpiexec  # requires PETSC_ARCH

# set CLAW_MPIFC to the proper Fortran compiler to use for MPI code
# e.g. mpif90 if that is defined in your shell, or gfortran *might* work.
# This will over-rule any FC environment variable.
export CLAW_MPIFC=mpif90
