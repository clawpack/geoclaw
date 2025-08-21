
.. _geoclaw_examples_bouss_radial_flat:

Radially symmetric problem on flat bottom with Boussinesq solver
=================================================================

A Gaussian hump of water at the origin spreads out over a flat bottom,
as specified by the topo file `flat100.txt` (uniform water depth 100 m).
The equations are solved in Cartesian coordinates with the SGN equations.

Running the GeoClaw Boussinesq solvers requires PETSc and MPI.
For more details see the documentation
  https://www.clawpack.org/bouss2d.html
and
  $CLAW/geoclaw/examples/bouss/README.txt
Run
  make check
in this directory to check if things seem ok for running this code.

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
