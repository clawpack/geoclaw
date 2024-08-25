
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools

try:
    TG32412 = np.loadtxt('32412_notide.txt')
except:
    print("*** Could not load DART data file")

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
    plotdata.format = 'ascii'    # 'ascii' or 'binary' to match setrun.py


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface at time h:m:s'
    plotaxes.aspect_latitude = -30.
    plotaxes.xticks_fontsize = 10
    plotaxes.yticks_fontsize = 10
    plotaxes.xlabel = 'longitude'
    plotaxes.ylabel = 'latitude'

    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.8
    plotitem.colorbar_extend = 'both'
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [-120,-60]
    plotaxes.ylimits = [-60,0]

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
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = 1/3600.  # convert to hours
    plotaxes.xlimits = [0, 9]
    plotaxes.ylimits = [-0.3, 0.3]
    plotaxes.title = 'Surface'
    plotaxes.title_fontsize = 15
    plotaxes.time_label = 'time (hours)'
    plotaxes.ylabel = 'surface elevation (m)'
    plotaxes.grid = True

    def add_obs(current_data):
        from pylab import plot, legend
        gaugeno = current_data.gaugeno
        if gaugeno == 32412:
            try:
                tgauge = TG32412[:,0] / 3600.  # convert to hours
                plot(tgauge, TG32412[:,1], 'r', linewidth=1.0)
                legend(['GeoClaw','Obs'],loc='lower right')
            except: pass

    plotaxes.afteraxes = add_obs

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3  # eta
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth':1.8}


    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False  # not being used

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'



    #-----------------------------------------
    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        from clawpack.visclaw import plot_timing_stats
        import os,sys
        try:
            timing_plotdir = plotdata.plotdir + '/_timing_figures'
            os.system('mkdir -p %s' % timing_plotdir)
            # adjust units for plots based on problem:
            units = {'comptime':'seconds', 'simtime':'hours', 
                     'cell':'millions'}
            plot_timing_stats.make_plots(outdir=plotdata.outdir, 
                                          make_pngs=True,
                                          plotdir=timing_plotdir, 
                                          units=units)
        except:
            print('*** Error making timing plots')

    otherfigure = plotdata.new_otherfigure(name='timing plots',
                    fname='_timing_figures/timing.html')
    otherfigure.makefig = make_timing_plots


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
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

