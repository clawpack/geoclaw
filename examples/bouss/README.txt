
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
or directly in the application Makefile, e.g. see the lines commented out in
radial_flat/Makefile.

A file petscMPIoptions is also required to set some PETSc parameters for the
iterative method used to solve the large sparse linear systems that arise at 
each refinement level when the Boussinesq equations are solved. 
One of the environment variables mentioned above points to this file, and a
sample is included in this directory.
