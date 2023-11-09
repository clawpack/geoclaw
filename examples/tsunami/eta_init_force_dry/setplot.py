
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools
import os,sys

cmax = 0.5
cmin = -cmax

cmax_land = 2.

#--------------------------
def setplot(plotdata=None):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'



    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr

    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('%s after earthquake' % timestr)


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Computational domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True


    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, savefig
        gca().set_aspect(1.)
        title_hours(current_data)
        ticklabel_format(useOffset=False)
        xticks(rotation=20)


    plotaxes.afteraxes = aa


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,0,0,0]
    plotitem.amr_data_show = [1,1,1,1,1,0,0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax_land
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,0,0]
    plotitem.amr_data_show = [1,1,1,1,1,0,0]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [1,0,0]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    #-----------------------------------------
    # Figure for coastal area
    #-----------------------------------------
    x1,x2,y1,y2 = [-0.005, 0.016, -0.01, 0.01]

    plotfigure = plotdata.new_plotfigure(name="coastal area", figno=11)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize': (6,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.scaled = False

    plotaxes.xlimits = [x1, x2]
    plotaxes.ylimits = [y1, y2]

    def aa_withbox(current_data):
        from pylab import plot
        x1,x2,y1,y2 = (-0.009259, 0.013796, -0.005093, 0.005000)
        if current_data.t > 5*60.:
            plot([x1,x1,x2,x2,x1], [y1,y2,y2,y1,y1], 'w--')
        aa(current_data)
        
    plotaxes.afteraxes = aa_withbox

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.4
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax_land
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show = 0

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-2,-1,0,1,2]
    plotitem.amr_contour_colors = ['yellow']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [0,0,0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0



    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        import os
        from clawpack.visclaw import plot_timing_stats
        try:
            timing_plotdir = plotdata.plotdir + '/_timing_figures'
            os.system('mkdir -p %s' % timing_plotdir)
            units = {'comptime':'minutes', 'simtime':'minutes', 'cell':'millions'}
            plot_timing_stats.make_plots(outdir=plotdata.outdir, make_pngs=True,
                                          plotdir=timing_plotdir, units=units)
            os.system('cp %s/timing.* %s' % (plotdata.outdir, timing_plotdir))
        except:
            print('*** Error making timing plots')

    otherfigure = plotdata.new_otherfigure(name='timing',
                    fname='_timing_figures/timing.html')
    otherfigure.makefig = make_timing_plots


    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'        # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata
