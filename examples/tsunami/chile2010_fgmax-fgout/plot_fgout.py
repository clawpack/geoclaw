
from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta 

fgno = 1
outdir = '_output'
output_format = 'binary'

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, output_format)
fgout_grid.read_fgout_grids_data()

# Plot one frame of fgout data 
fgframe = 20
fgout = fgout_grid.read_frame(2)

fgout = fgout_grid.read_frame(fgframe)

figure(1, figsize=(8,8))
imshow(flipud(fgout.B.T), extent=fgout.extent_edges,
       cmap=geoplot.land_colors)

clim(0,100)

eta_water = where(fgout.h > 0, fgout.eta, nan)
imshow(flipud(eta_water.T), extent=fgout.extent_edges,
       cmap=geoplot.tsunami_colormap)

clim(-0.2, 0.2)

title('Surface at time %s' % timedelta(seconds=fgout.t))

fname = 'fgout_frame%s.png' % str(fgframe).zfill(4)
savefig(fname)
print('Created ',fname)

