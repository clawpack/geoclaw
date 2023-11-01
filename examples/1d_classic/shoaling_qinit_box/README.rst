
.. _geoclaw_1d/examples/shoaling_qinit_box:

Shoaling example with square wave data
========================================

The initial wave is specified in `qinit.f90` as a Gaussian perturbation in
the sea surface `eta`, with the momentum set to `eta * sqrt(g*h)` so that
the wave is initially purely right-going on the constant-depth ocean.

This example was adapted from the example at
https://github.com/rjleveque/shoaling_paper_figures/tree/master/qinit_box,
originally developed to generate figures 2 and 3 for the paper::

    Shoaling on Steep Continental Slopes: Relating Transmission and 
    Reflection Coefficients to Green's Law, 
    by J. George, D. I. Ketcheson, and R. J. LeVeque,
    Pure and Applied Geophysics, 2019.
    DOI 10.1007/s00024-019-02316-y

See http://faculty.washington.edu/rjl/pubs/Shoaling2019 for additional links.

The file `make_celledges.py` sets up the domain and computational grid.
A piecewise linear topography is defined by specifying the topography `z`
value at a set of nodes `x` in the `xzpairs` list.  Set up to have an
ocean of depth 3200m and a shelf of depth 200m.

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


Version
-------

Updated when merged into geoclaw, November 2023

