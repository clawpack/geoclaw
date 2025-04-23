"""
Make an mp4 animation of fgout grid results that also includes a transect
of the solution, as often desired in animations.

This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.

Make the animation via:
    python make_fgout_animation_transect.py

If this script is executed in IPython or a notebook it may go into
an infinite loop for reasons unknown.  If so, close the figure to halt.

To view individual fgout frames interactively, this should work:
    import make_fgout_animation_transect
    make_fgout_animation_transect.update(fgframeno)  # for desired fgout frameno

Uses blit==False so that update_artists tuple does not need to be returned
from the update function.
"""

import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
import os, glob
from clawpack.visclaw import plottools, geoplot, gridtools, animation_tools
from clawpack.geoclaw import util
from matplotlib import animation, colors
from datetime import timedelta

from clawpack.geoclaw import fgout_tools
    
fgno = 1  # which fgout grid

outdir = '_output'

if 1:
    # use all fgout frames in outdir:
    fgout_frames = glob.glob(os.path.join(outdir, \
                                          'fgout%s.t*' % str(fgno).zfill(4)))
    nout = len(fgout_frames)
    fgframes = range(1, nout+1)
    print('Found %i fgout frames in %s' % (nout,outdir))
else:
    # set explicitly, e.g. to test with only a few frames
    fgframes = range(1,26)  # frames of fgout solution to use in animation

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir) 
fgout_grid.read_fgout_grids_data()


# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

fgout = fgout_grid.read_frame(fgframes[0])

plot_extent = fgout.extent_edges
ylat = fgout.Y.mean()  # for aspect ratio of plots

fig = figure(figsize=(12,7))

# ---------------------------------
# axis for planview plot of ocean:
ax = axes([.1,.1,.4,.8])

ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])

plottools.pcolorcells(fgout.X,fgout.Y,fgout.B, cmap=geoplot.land_colors)
clim(0,100)

eta = ma.masked_where(fgout.h<0.001, fgout.eta)

eta_plot = plottools.pcolorcells(fgout.X,fgout.Y,eta,
                                 cmap=geoplot.tsunami_colormap)
clim(-0.3,0.3)
cb = colorbar(eta_plot, extend='both', shrink=0.5,
              orientation='horizontal', anchor=(0.4,1))
cb.set_label('meters')
title_text = title('Surface at time %s' % timedelta(seconds=fgout.t))

ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])

# ---------------------------------
# Transect:

x1trans = -110.; y1trans = -20.
x2trans = -72.;  y2trans = -35.

plot([x1trans,x2trans], [y1trans,y2trans],'k-',linewidth=0.8)
text(x1trans+1,y1trans+1,'Transect', ha='center', fontsize=10)

# define points on transect:
npts = 1000  # number of points on transect
if 0:
    # straight line on longitude-latitude plane:
    xtrans = linspace(x1trans,x2trans,npts)
    ytrans = linspace(y1trans,y2trans,npts)
else:
    # great circle on earth:
    xtrans,ytrans = util.gctransect(x1trans,y1trans,x2trans,y2trans,npts)
    

def extract_transect(fgout_soln,xtrans,ytrans):
    """
    Interpolate from fgout_solution to points on the transect, taking value
    from cell the fgout point lies in, giving piecewise constant interpolant
    when npts is large.
    Alternatively could use method='linear' for linear interpolation.
    """

    eta1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                   fgout_soln.eta, xtrans, ytrans,
                                   method='nearest')
    B1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                 fgout_soln.B, xtrans, ytrans,
                                 method='nearest')
    eta1d = where(B1d>0, nan, eta1d)   # mask onshore region
    return B1d, eta1d

# ----------------------------------------
# Axes for transect of surface elevation:

axtrans1 = axes([.55,.6,.4,.3])

axtrans1.set_title('Surface on transect')

axtrans1.set_xlim(x1trans,x2trans)
axtrans1.set_ylim(-1,1)

Btrans1, etatrans1 = extract_transect(fgout,xtrans,ytrans)
#import pdb; pdb.set_trace()

Btrans, etatrans = Btrans1, etatrans1

# surface plot:
etatrans1_plot, = axtrans1.plot(xtrans, etatrans, 'b')
axtrans1.grid(True)

# ----------------------------------------
# Axes for transect of topography:

axtrans2 = axes([.55,.1,.4,.3])

axtrans2.set_title('Topography on transect')

axtrans2.set_xlim(x1trans,x2trans)
axtrans2.set_ylim(-5000,1000)

Btrans, etatrans = extract_transect(fgout,xtrans,ytrans)

# filled regions:
Bfill_plot = axtrans2.fill_between(xtrans, Btrans-1e5, Btrans, 
                                   color=[.7,1,.7,1]) # light green
etafill_plot = axtrans2.fill_between(xtrans, Btrans, etatrans, 
                                     color=[.7,.7,1,1])  # light blue

# surface and topo solid lines:
etatrans2_plot, = axtrans2.plot(xtrans, etatrans, 'b')
Btrans2_plot, = axtrans2.plot(xtrans, Btrans, 'g')


# create a dummy figure and axes, only needed to update fill_between artists:
figdummy,axdummy = subplots()


def update(fgframeno):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    Assumes blit==False in call to animation.FuncAnimation below,
    so tuple of update_artists does not need to be returned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))
        
    # reset title to current time:
    title_text.set_text('Surface at time %s' % timedelta(seconds=fgout.t))

    # reset surface eta to current state:
    eta = ma.masked_where(fgout.h<0.001, fgout.eta)
    eta_plot.set_array(eta.T.flatten())
        
    # update transect data:
    Btrans, etatrans = extract_transect(fgout,xtrans,ytrans)    
    
    # update lines plotted:
    etatrans1_plot.set_data(xtrans,etatrans)
    etatrans2_plot.set_data(xtrans,etatrans)
    Btrans2_plot.set_data(xtrans,Btrans)

    # update the PolyCollections for fill_between plots:
    # There doesn't seem to be an easier way to do this...             
    dummy = axdummy.fill_between(xtrans, Btrans-1e5, Btrans, color=[.5,1,.5,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    Bfill_plot.set_paths([dp.vertices])

    dummy = axdummy.fill_between(xtrans, Btrans, etatrans, color=[.5,.5,1,1])
    dp = dummy.get_paths()[0]
    dummy.remove()
    etafill_plot.set_paths([dp.vertices])
    

if __name__ == '__main__':

    print('Making anim...')
    anim = animation.FuncAnimation(fig, update, frames=fgframes, 
                                   interval=200, blit=False)
    
    # Output files:
    name = 'fgout_animation_with_transect'

    fname_mp4 = name + '.mp4'

    #fname_html = None
    fname_html = name + '.html'

    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        animation_tools.make_html(anim, file_name=fname_html, title=name)
