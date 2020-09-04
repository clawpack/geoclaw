#!/usr/bin/env python
# coding: utf-8

# # process_fgmax
# 
# Read in fgmax results and produce plots.

# In[ ]:


from pylab import *


# In[ ]:


import os,sys
import glob
from scipy.interpolate import RegularGridInterpolator
import matplotlib as mpl
from matplotlib import colors

from clawpack.geoclaw import topotools, dtopotools
from clawpack.visclaw import colormaps
from clawpack.visclaw.plottools import pcolorcells
from clawpack.geoclaw import fgmax_tools


# In[ ]:


save_figs = True
fgmax_plotdir = '_plots/fgmax_plots'


# In[ ]:


os.system('mkdir -p %s' % fgmax_plotdir)
def savefigp(fname):
    global save_figs
    if save_figs:
        fullname = '%s/%s' % (fgmax_plotdir, fname)
        savefig(fullname)
        print('Created ', fullname)
    else:
        print('save_figs = False')


# In[ ]:


outdir = '_output'
t_files = glob.glob(outdir + '/fort.t0*')
times = []
for f in t_files:
    lines = open(f,'r').readlines()
    for line in lines:
        if 'time' in line: 
            t = float(line.split()[0])
    times.append(t)
times.sort()
print('Output times found: ',times)
if len(times) > 0:
    t_hours = times[-1] / 3600.
    print('\nfgmax results are presumably from final time: %.1f seconds = %.2f hours'          % (times[-1], t_hours))
else:
    t_hours = nan


# In[ ]:


# Read fgmax data:
fgno = 1
fg = fgmax_tools.FGmaxGrid()
fg.read_fgmax_grids_data(fgno)

fg.read_output(outdir=outdir)


# In[ ]:


zmin = -60.
zmax = 20.
land_cmap = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                     0.25:[0.0,1.0,0.0],
                                      0.5:[0.8,1.0,0.5],
                                      1.0:[0.8,0.5,0.2]})

sea_cmap = colormaps.make_colormap({ 0.0:[0,0,1], 1.:[.8,.8,1]})

cmap, norm = colormaps.add_colormaps((land_cmap, sea_cmap),
                                     data_limits=(zmin,zmax),
                                     data_break=0.)                                   



figure(figsize=(8,8))
pc = pcolorcells(fg.X, fg.Y, fg.B, cmap=cmap, norm=norm)  

cb = colorbar(pc,shrink=0.5,extend='both')
cb.set_label('meters')
cb.set_ticks(hstack((linspace(zmin,0,5), linspace(0,zmax,5))))

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20);
    
title('GeoClaw B topography on fg1 grid');


# In[ ]:


fg.B0 = fg.B  # no seafloor deformation in this problem
fg.h_onshore = ma.masked_where(fg.B0 < 0., fg.h)


# In[ ]:


bounds_depth = array([1e-6,0.5,1.0,1.5,2,2.5,3.0])


cmap_depth = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],                 [1,.7,.7], [1,.4,.4], [1,0,0]])

# Set color for value exceeding top of range to purple:
cmap_depth.set_over(color=[1,0,1])

# Set color for land points without inundation to light green:
cmap_depth.set_under(color=[.7,1,.7])

norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)
    

figure(figsize=(8,8))
pc = pcolorcells(fg.X, fg.Y, fg.h_onshore, cmap=cmap_depth, norm=norm_depth)
cb = colorbar(pc, extend='max', shrink=0.7)
cb.set_label('meters')
contour(fg.X, fg.Y, fg.B, [0], colors='g')

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
title('Maximum Onshore flow depth over %.2f hours\nfgmax grid %s' % (t_hours,fgno))
savefigp('fgmax%s_h_onshore.png' % str(fgno).zfill(4))


# In[ ]:


bounds_speed = np.array([1e-6,0.5,1.0,1.5,2,2.5,3,4.5,6])
cmap_speed = mpl.colors.ListedColormap([[.9,.9,1],[.6,.6,1],                 [.3,.3,1],[0,0,1], [1,.8,.8],                 [1,.6,.6], [1,.3,.3], [1,0,0]])


bounds_speed = np.array([1e-6,0.5,1.0,1.5,2,2.5,3,4.5])
cmap_speed = mpl.colors.ListedColormap([[.9,.9,1],[.6,.6,1],                 [.3,.3,1],[0,0,1], [1,.8,.8],                 [1,.6,.6], [1,0,0]])

# Set color for value exceeding top of range to purple:
cmap_speed.set_over(color=[1,0,1])

# Set color for land points without inundation to light green:
cmap_speed.set_under(color=[.7,1,.7])

norm_speed = colors.BoundaryNorm(bounds_speed, cmap_speed.N)

figure(figsize=(8,8))
pc = pcolorcells(fg.X, fg.Y, fg.s, cmap=cmap_speed, norm=norm_speed)
cb = colorbar(pc, extend='max', shrink=0.7)
cb.set_label('m/s')
contour(fg.X, fg.Y, fg.B0, [0], colors='g')

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
title('Maximum speed over %.2f hours\nfgmax grid %s' % (t_hours,fgno))
savefigp('fgmax%s_speed.png' % str(fgno).zfill(4))


# Save this so we can plot the topo below...

# In[ ]:


import copy
fg1 = copy.copy(fg)


# ## Read fgmax values specified on a Transect

# In[ ]:


# Read fgmax data:
fgno = 2
fg = fgmax_tools.FGmaxGrid()
fg.read_fgmax_grids_data(fgno)
fg.read_output(outdir=outdir)
xx = fg.X
yy = fg.Y


# In[ ]:


figure(figsize=(8,8))
pc = pcolorcells(fg1.X, fg1.Y, fg1.B, cmap=cmap, norm=norm)  

cb = colorbar(pc,shrink=0.5,extend='both')
cb.set_label('meters')
cb.set_ticks(hstack((linspace(zmin,0,5), linspace(0,zmax,5))))

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20);
    
plot(xx,yy,'r')
title('GeoClaw B topography values on fg1 grid\n with transect from fg2');


# In[ ]:


figure(figsize=(12,4))
fill_between(xx, fg.B, fg.B+fg.h, color=[.5,.5,1])
plot(xx,fg.B+fg.h,'b')
plot(xx,fg.B,'g')
plot(xx, ma.masked_where(fg.B>0, 0*xx), 'k')
grid(True)
ylim(-10,20);
title('Maximum elevation over %.2f hours\nfgmax grid %s' % (t_hours,fgno))
savefigp('fgmax%s_surface.png' % str(fgno).zfill(4));


# ## Read fgmax points as specified on a masked grid

# In[ ]:


fgno = 3
fg = fgmax_tools.FGmaxGrid()
fg.read_fgmax_grids_data(fgno)

fg.read_output(outdir=outdir)


# In[ ]:


fg.B0 = fg.B  # no seafloor deformation in this problem
fg.h_onshore = ma.masked_where(fg.B0 < 0., fg.h)


# In[ ]:


figure(figsize=(8,8))
pc = pcolorcells(fg.X, fg.Y, fg.B, cmap=cmap, norm=norm)
cb = colorbar(pc, extend='both', shrink=0.7)
cb.set_label('meters')
cb.set_ticks(hstack((linspace(zmin,0,5), linspace(0,zmax,5))))

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
title('GeoClaw B at points selected as fgmax grid\nfgmax grid %s' % fgno);


# In[ ]:


figure(figsize=(8,8))
pc = pcolorcells(fg.X, fg.Y, fg.h_onshore, cmap=cmap_depth, norm=norm_depth)
cb = colorbar(pc, extend='max', shrink=0.7)
cb.set_label('meters')
contour(fg.X, fg.Y, fg.B0, [0], colors='g')

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
title('Maximum Onshore flow depth over %.2f hours' % t_hours);
savefigp('fgmax%s_h_onshore.png' % str(fgno).zfill(4))


# In[ ]:


figure(figsize=(8,8))
pc = pcolorcells(fg.X, fg.Y, fg.s, cmap=cmap_speed, norm=norm_speed)
cb = colorbar(pc, extend='max', shrink=0.7)
cb.set_label('m/s')
contour(fg.X, fg.Y, fg.B0, [0], colors='g')

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
title('Maximum speed over %.2f hours\nfgmax grid %s' % (t_hours,fgno))
savefigp('fgmax%s_speed.png' % str(fgno).zfill(4))


# ### View fgmax points selected
# 
# This isn't generally needed, but if you want to inspect the file that specified fgmax points originally:

# In[ ]:


fg3input = topotools.Topography(path=fg.xy_fname, topo_type=3)
fg3input.X.shape

figure(figsize=(8,8))
pc = pcolorcells(fg3input.X, fg3input.Y, fg3input.Z)
cb = colorbar(pc, shrink=0.7)

gca().set_aspect(1./cos(48*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20);


# ## Read points with `point_style == 0`

# In[ ]:


# Read fgmax data:
fg = fgmax_tools.FGmaxGrid()
fg.read_fgmax_grids_data(4)

fg.read_output(outdir=outdir)
print('\n      x          y       max depth')
for j in range(fg.npts):
    print('%10.3f %10.3f %10.3f'  % (fg.X[j], fg.Y[j], fg.h[j]))


# In[ ]:


# Read fgmax data:
fg = fgmax_tools.FGmaxGrid()
fg.read_fgmax_grids_data(5)

fg.read_output(outdir=outdir)
print('\n      x          y       max speed')
for j in range(fg.npts):
    print('%10.3f %10.3f %10.3f'  % (fg.X[j], fg.Y[j], fg.s[j]))


# In[ ]:




