
.. _geoclaw_examples_tsunami_island-particles:

Flow past an island with particle tracking
==========================================

Lagrangian gauges (introduced in v5.7.0) are used to track passive particles 
carried by the flow.  In this example an advancing hydraulic jump is used to
set up flow past a conical island.  Particle tracking helps to visualize the
vortices generated behind the island.

To create the topo file before running the code::

    make topo


In this code, :math:`x` and :math:`y` are in meters (coordinate_system=1 
in `setrun.py`).

In `setrun.py` 100 gauges are specified, 20 stationary and the other 
80 are lagrangian.

In `setplot.py`, the `clawpack.visclaw.particle_tools` module is used to
plot particle locations and particle paths.

