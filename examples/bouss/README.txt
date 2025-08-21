
GeoClaw Boussinesq solver examples (in 2D)

Some 1D Boussinesq solver examples are in $CLAW/geoclaw/examples/1d_classic
in directories starting with bouss_

The radial_flat subdirectory contains one example using the Boussinesq
solver introduced in Clawpack v5.10.0.

Several other examples can be found in this repository:
    https://github.com/rjleveque/ImplicitAMR-paper
which are described in the paper found at:
    Implicit Adaptive Mesh Refinement for Dispersive Tsunami Propagation
    by Marsha J. Berger and Randall J. LeVeque, submitted, 2023
    Preprint: https://arxiv.org/abs/2307.05816

Running these codes requires PETSc  Version 3.20 (or later) in order to use
OpenMP along with MPI.  Some flags have to be set as environment variables
or directly in the application Makefile.

For more information on installing PETSc and setting these environment
variables properly, see:
  https://www.clawpack.org/bouss2d.html


**Update:** Clawpack 5.12.0 now requires PETSc Version 3.23 (or later).

The file setenv.sh illustrates how you might set some environment
variables for the bash shell.

A file petscMPIoptions is also required to set some PETSc parameters for the
iterative method used to solve the large sparse linear systems that arise at 
each refinement level when the Boussinesq equations are solved. 
The environment variable PETSC_OPTIONS should point to the version of this
file that you wish to use, see setenv.sh for an example of setting this to
point to the sample petscMPIoptions file included in this directory.


**Check:**

The example Makefile in
    $CLAW/geoclaw/examples/bouss/radial_flat/Makefile
contains an option so that you can do 
    $ make check
at the command line to see how various options are set and to check that
they look correct (or help debug if the code does not run).
This is an enhancement of the 'make check' option added to the common
Makefile for all Clawpack applications in v5.12.0, and should produce
something like this:

===================
CLAW = /Users/rjl/git/clawpack
OMP_NUM_THREADS = 6
BOUSS_MPI_PROCS = 6
PETSC_OPTIONS=-options_file /Users/rjl/git/clawpack/geoclaw/examples/bouss/petscMPIoptions
PETSC_DIR = /Users/rjl/git/Clones/petsc
PETSC_ARCH = arch-darwin-c-opt
CLAW_MPIEXEC = mpiexec
RUNEXE = mpiexec -n 6
EXE = /Users/rjl/git/clawpack/geoclaw/examples/bouss/radial_flat/xgeoclaw
CLAW_MPIFC = mpif90
FC = mpif90
FFLAGS = -O2 -fopenmp -DHAVE_PETSC -ffree-line-length-none
LFLAGS = -L/Users/rjl/git/Clones/petsc/arch-darwin-c-opt/lib -lpetsc -fopenmp
OUTDIR = _output
PLOTDIR = _plots
===================

