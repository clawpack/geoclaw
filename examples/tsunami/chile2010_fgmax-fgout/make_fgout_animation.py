"""
Make an mp4 animation of fgout grid results. 
This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.

Make the animation via:
    python make_fgout_animation.py

If this script is executed in IPython or a notebook it may go into
an infinite loop for reasons unknown.  If so, close the figure to halt.

To view individual fgout frames interactively, this should work:
    import make_fgout_animation
    make_fgout_animation.update(fgframeno)  # for desired fgout frame no

"""

import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
import os, glob
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
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

fgout1 = fgout_grid.read_frame(fgframes[0])

plot_extent = fgout1.extent_edges
ylat = fgout1.Y.mean()  # for aspect ratio of plots

fig,ax = subplots(figsize=(8,7))

ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])

plottools.pcolorcells(fgout1.X,fgout1.Y,fgout1.B, cmap=geoplot.land_colors)
clim(0,100)

eta = ma.masked_where(fgout1.h<0.001, fgout1.eta)

eta_plot = plottools.pcolorcells(fgout1.X,fgout1.Y,eta,
                                 cmap=geoplot.tsunami_colormap)
clim(-0.3,0.3)
cb = colorbar(eta_plot, extend='both', shrink=0.7)
cb.set_label('meters')
title_text = title('Surface at time %s' % timedelta(seconds=fgout1.t))

ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])


blit = False
if blit:
    # The artists that will be updated for subsequent frames:
    update_artists = (eta_plot, title_text)
    
# Note that update_artists is only needed if blit==True in call to
# animation.FuncAnimation below.
# Using blit==False does not seem to slow down creation of the animation by much
# and slightly simplifies modification of this script to situations where more
# artists are updated.
        
def update(fgframeno):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    Note: Even if blit==True in call to animation.FuncAnimation,
    the update_artists do not need to be passed in, unpacked, and repacked
    as in an earlier version of this example (Clawpack version <= 5.10.0).
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))
        
    # reset title to current time:
    title_text.set_text('Surface at time %s' % timedelta(seconds=fgout.t))

    # reset surface eta to current state:
    eta = ma.masked_where(fgout.h<0.001, fgout.eta)
    eta_plot.set_array(eta.T.flatten())
        
    if blit:
        return update_artists


if __name__ == '__main__':

    print('Making anim...')
    anim = animation.FuncAnimation(fig, update, frames=fgframes, 
                                   interval=200, blit=blit)
    
    # Output files:
    name = 'fgout_animation'

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

