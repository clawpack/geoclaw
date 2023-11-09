"""
Make an mp4 animation of fgout grid results. 
This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.
The tuple update_artists contains the list of Artists that must be changed
in update.  Modify this as needed.

"""

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from datetime import timedelta

from clawpack.geoclaw import fgout_tools
    
fgno = 1  # which fgout grid

outdir = '_output'
format = 'binary'  # format of fgout grid output

fgframes = range(1,26)  # frames of fgout solution to use in animation

figsize = (10,8)

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 


# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

fgout1 = fgout_grid.read_frame(fgframes[0])

plot_extent = fgout1.extent_edges

ylat = fgout1.Y.mean()  # for aspect ratio of plots

fig,ax = subplots(figsize=figsize)

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


# The artists that will be updated for subsequent frames:
update_artists = (eta_plot, title_text)
        
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))
    
    # unpack update_artists (must agree with definition above):
    eta_plot, title_text = update_artists
        
    # reset title to current time:
    title_text.set_text('Surface at time %s' % timedelta(seconds=fgout.t))

    # reset surface eta to current state:
    eta = ma.masked_where(fgout.h<0.001, fgout.eta)
    eta_plot.set_array(eta.T.flatten())
        
    update_artists = (eta_plot, title_text)
    return update_artists

def plot_fgframe(fgframeno):
    """
    Convenience function for plotting one frame.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=fgframes, 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':

    anim = make_anim()
    
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
    
    
    
    
