
.. _geoclaw_examples_bouss_radial_flat:

Radially symmetric problem on flat bottom with Boussinesq solver
=================================================================

A Gaussian hump of water at the origin spreads out over a flat bottom,
as specified by the topo file `flat100.txt` (uniform water depth 100 m).
The equations are solved in Cartesian coordinates with the SGN equations.

Non-PETSc version:
------------------

This version is new in v5.15.0 and uses a stand-alone GMRES implementation,
which does not require PETSc or MPI.

This is now the default for this example and can be compiled with the
Makefile, which is similar to all other standard GeoClaw Makefiles.

Set 
    rundata.bouss_data.bouss_solver = 1 
in setrun.py to use this version.


PETSc version:
--------------

This is the version from v5.14.0 and earlier, and requires PETSc and MPI.

To use this version, compile and run with Makefile_petsc, e.g.
    make new -f Makefile_petsc
    make .output -f Makefile_petsc

Set 
    rundata.bouss_data.bouss_solver = 3 
in setrun.py to use this version.

For more details see the documentation
  https://www.clawpack.org/bouss2d.html
and
  $CLAW/geoclaw/examples/bouss/README.txt
Run
  make check
in this directory to check if things seem ok for running this code.


Notes on the problem setup:
---------------------------

A flagregion specified by `RuledRectangle_Diagonal.data` (that is created by
code in `setrun.py` is used to allow flagging for refinement to
level 2 only near the diagonal (for abs(x-y) < 1000).  The code is set up to
run with 2 levels and captures the waves well except for the highest
frequency trailing wave.  

Adding a third level of refinement near the origin improves the accuracy of
the trailing wave. This can be done by setting ::

    amrdata.amr_levels_max = 3

in `setrun.py` and then the flagregion near the origin comes into play.

Before plotting, run the code in the `1d_radial` subdirectory to
create a reference solution on a fine grid (5000 cells).



Version
-------

- Created for v5.10.0
- Modified for v5.15.0
