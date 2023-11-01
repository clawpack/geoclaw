
.. _geoclaw_1d/examples/bouss_wavetank_matsuyama:

Wave tank simulation compared to Matsuyama et al. (2007) experiment
===================================================================

This example uses the wave tank dimensions specified for Case 024 in ::

    A Study of Tsunami Wave Fission in an Undistorted Experiment,
    by M. Matsuyama, M. Ikeno, T. Sakiyama, and T. Takeda, 
    Pure Appl. Geophys. 164 (2007), pp. 617-631.
    DOI 10.1007/s00024-006-0177-0

Boundary conditions are implemented in `bc1.f` that mimic the wave maker.

This example is set up to use the Boussinesq equations, and reproduces the
break up of the dispersive wave into a train of solitary waves.

The file `make_celledges.py` sets up the domain and computational grid.
A piecewise linear topography is defined by specifying the topography `z`
value at a set of nodes `x` in the `xzpairs` list. 

A nonuniform grid with `mx` grid cells is used with cell widths related
to the still water depth in such a way that the Courant number is roughly
constant in deep water and onto the shelf, and with uniform grid cells
near shore and onshore where the water depth is less than `hmin`.

Executing `make_celledges.py` creates a file `celledges.data` that contains
the cell edges.  This file must be created before running GeoClaw.

In GeoClaw a mapped grid is used with a `mapc2p` function specified in
`setrun.py` that is generated from the `celledges.data`.  The computational
grid specified in `setrun.py` is always `0 <= xc <= 1`.  Set::

    rundata.grid_data.grid_type = 2
    
to indicate a mapped grid.

In this example the physical `x` coordiate is in meters, set by specifying::

    rundata.geo_data.coordinate_system = 1

To use::

    make topo     # executes make_celledges.py
    make .output  # compile, make data, and run
    make .plots   # to create _plots (or plot interactively with Iplotclaw)

    python plot_gauges.py  # to create a plot of gauges to compare to paper

The `plot_gauges.py` script should create a plot similar to
`GeoClawFigure5.png <GeoClawFigure5.png>`__ 
that can be compared to 
`Figure 5 <MatsuyamaFigure5.png>`__ 
in the paper.

The gauge plot produced also includes wave tank observations (as a red
curve) for some gauges.  This data comes from the file
data_wavegauge.csv (kindly provided by Prof. Matsuyama).

Version
-------

Updated when merged into geoclaw, November 2023


