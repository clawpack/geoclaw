"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw import nonuniform_grid_tools

xlower = -160.8
xupper = 40.

xzpairs = [(-160.8, -4),       # left edge
           (-125.5, -4),       # start of first slope
           (-90,    -0.45),    # start of beach
           (  0,     0),       # shore
           ( 40,     2)]       # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 8000
hmin = 0.05  # use uniform grid in shallower water
#hmin = 5  # try forcing a uniform grid everywhere
#hmin = 1  # try uniform grid near start of shelf

nonuniform_grid_tools.make_celledges_cfl(xlower, xupper, mx, topo_fcn,
        hmin, fname='celledges.data', plot_topo=True)

