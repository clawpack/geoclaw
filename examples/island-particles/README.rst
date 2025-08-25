
.. _geoclaw_examples_tsunami_island-particles:

Flow past an island with particle tracking
==========================================

Lagrangian gauges (introduced in v5.7.0) are used to track passive particles 
carried by the flow.  In this example an advancing hydraulic jump is used to
set up flow past a conical island.  Particle tracking helps to visualize the
vortices generated behind the island.


To create the topo file before running the code::

    make topo

This creates `island.tt3` (with topo_type == 3).  It also creates 
`qinit.xyz`, containing initial data for a surface perturbation
corresponding to a "dam break" problem at `x = 10`, leading to flow in the
positive `x` direction.


In this code, :math:`x` and :math:`y` are in meters (coordinate_system=1 
in `setrun.py`).

In `setrun.py` 100 gauges are specified, 20 stationary and the other 
80 are lagrangian.

Note that a very large `regrid_interval` is specified so that the refinement
patches are frozen and do not adapt to the flow.   This is done to
illustrate that particles are properly tracked when they move between grid
patches at different resolutions.  Three levels are used with the finest
level only around the island, as specified in the `regions` set up in
`setrun.py`.

In `setplot.py`, the `clawpack.visclaw.particle_tools` module is used to
plot particle locations and particle paths.

