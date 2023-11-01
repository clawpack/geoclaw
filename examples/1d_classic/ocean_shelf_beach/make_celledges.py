"""
Set up the domain and computational grid.
A piecewise linear topography is defined by specifying the topography `z`
value at a set of nodes `x` in the `xzpairs` list.

A nonuniform grid with `mx` grid cells is used with cell widths related
to the still water depth in such a way that the Courant number is roughly
constant in deep water and onto the shelf, and with uniform grid cells
near shore and onshore where the water depth is less than `hmin`.
"""

from pylab import *
from clawpack.geoclaw import nonuniform_grid_tools

x1 = -100e3
x2 = 1e3

xzpairs = [(-100e3,-4000),   # left edge
           ( -70e3,-4000),   # start of continental slope
           ( -20e3,-200),    # start of continental shelf
           (  -2e3,-200),    # start of beach
           (     0,   0),       # shoreline
           (   1e3, 100)]     # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 5000   # number of grid cells
hmin = 20.  # mininum depth for varying cell widths

nonuniform_grid_tools.make_celledges_cfl(x1, x2, mx, topo_fcn,
        hmin=hmin, fname='celledges.data', plot_topo=True)

