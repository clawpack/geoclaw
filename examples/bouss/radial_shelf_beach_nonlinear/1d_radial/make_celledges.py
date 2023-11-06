"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw import nonuniform_grid_tools


mx = 10000  # desired number of grid cells 10000 => 5m on shore
print('mx = %i' % mx)

x0 = 0.
x0_shore = 120e3     # initial shoreline
x0_beach = 100e3       # start of beach
x0_shelf = 80e3     # start of continental shelf
x0_slope = 40e3      # start of continental slope
x1 = x0_shore + 4e3 

z0_ocean = -3000.     # depth of ocean = depth at x0_slope
z0_shelf = -100.      # depth at x0_shelf
z0_beach = -100.       # depth at x0_beach
z0_shore = 0.      # depth at x0_shore
slope_of_beach = (z0_beach - z0_shore) / (x0_beach - x0_shore)
z1 = z0_shore + (x1-x0_shore)*slope_of_beach

xzpairs = [(x0, z0_ocean),
           (x0_slope, z0_ocean),
           (x0_shelf, z0_shelf),
           (x0_beach, z0_shelf),
           (x0_shore, z0_shore),
           (x1, z1)]
           
topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

hmin = 50.  # use uniform grid in shallower water

xp,z = nonuniform_grid_tools.make_celledges_cfl(x0, x1, mx, topo_fcn,
            hmin, fname='celledges.data', plot_topo=True)

figure(98, figsize=(9,3))
clf()
fill_between(xp,where(z<0,z,nan),0.,color=[.5,.5,1])
plot(xp,z,'g')
grid(True)
title('Radial ocean topography',fontsize=15)
#xlabel('Distance from center (km)',fontsize=12)
ylabel('meters',fontsize=12)
xticks([0,40e3,80e3,100e3,120e3,126e3],[0,40,80,100,120,'km'],fontsize=12)
yticks(fontsize=10)
text(20e3,-3100,'deep ocean',va='top',ha='center',fontsize=12)
text(60e3,-3100,'continental slope',va='top',ha='center',fontsize=12)
text(90e3,-3100,'shelf',va='top',ha='center',fontsize=12)
text(110e3,-3100,'beach',va='top',ha='center',fontsize=12)
xlim(0,126e3)
ylim(-3400,400)

fname = 'radial_ocean_topo.pdf'
savefig(fname, bbox_inches='tight')
print('Created ',fname)
