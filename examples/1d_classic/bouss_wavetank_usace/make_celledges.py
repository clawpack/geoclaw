"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw import nonuniform_grid_tools

xlower = -15. #-0.98
xupper = 8.19

xzpairs = [(-0.98, -0.218),       # left edge
           (    0, -0.218),       # start of first slope
           ( 4.36, -0.1357),      # start of second slope
           ( 7.29, -0.1162),      # start of third slope
           ( 8.19, -0.0470)]      # right edge

# flat:
#xzpairs = [(xlower, -0.218),       # left edge
#           (xupper, -0.218)]       # right edge

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 800
hmin = 0.04  # use uniform grid in shallower water

nonuniform_grid_tools.make_celledges_cfl(xlower, xupper, mx, topo_fcn,
        hmin, fname='celledges.data', plot_topo=True)

