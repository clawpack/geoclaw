
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 


from __future__ import absolute_import
from __future__ import print_function

from clawpack.visclaw import gaugetools

from clawpack.visclaw import particle_tools
from clawpack.visclaw import legend_tools
from importlib import reload
reload(particle_tools)


#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotdata.format = 'ascii'                # Format of output

    print('Reading all gauges...')
    gauge_solutions = particle_tools.read_gauges(gaugenos='all', 
                                                 outdir=plotdata.outdir)

    gaugenos_lagrangian = [k for k in gauge_solutions.keys() \
                if gauge_solutions[k].gtype=='lagrangian']
    gaugenos_stationary = [k for k in gauge_solutions.keys() \
                if gauge_solutions[k].gtype=='stationary']

    print('+++ gaugenos_lagrangian: ',gaugenos_lagrangian)
    
    def add_particles(current_data):
        t = current_data.t

        # plot recent path:
        t_path_length = 10.   # length of path trailing particle
        kwargs_plot_path = {'linewidth':1, 'color':'k'}
        particle_tools.plot_paths(gauge_solutions, 
                                  t1=t-t_path_length, t2=t, 
                                  gaugenos=gaugenos_lagrangian, 
                                  kwargs_plot=kwargs_plot_path)

        # plot current location:
        kwargs_plot_point = {'marker':'o','markersize':3,'color':'k'}
        particle_tools.plot_particles(gauge_solutions, t, 
                                      gaugenos=gaugenos_lagrangian, 
                                      kwargs_plot=kwargs_plot_point)  

        # plot any stationary gauges:
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos_stationary, format_string='kx', add_labels=False)
        kwargs={'loc':'upper left'}
        legend_tools.add_legend(['Lagrangian particle','Stationary gauge'],
                linestyles=['',''], markers=['o','x'],
                loc='lower right', framealpha=0.5, fontsize=10)


    def speed(current_data):
        from pylab import sqrt, where, zeros
        from numpy.ma import masked_where, allequal
        q = current_data.q
        h = q[0,:,:]
        hs = sqrt(q[1,:,:]**2 + q[2,:,:]**2)
        s = where(h>1e-3, hs/h, 0.)
        s = masked_where(h<1e-3, s)
        s = s * 1.94384  # convert to knots
        return s

    speed_cmap = colormaps.make_colormap({0:[0,1,1], 0.5:[1,1,0], 1:[1,0,0]})


    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    plotfigure.kwargs = {'figsize': (9,4)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Speed'
    plotaxes.scaled = True
    plotaxes.xlimits = [0,80]
    plotaxes.ylimits = [0,50]
    plotaxes.afteraxes = add_particles

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = geoplot.surface_or_depth
    plotitem.plot_var = speed
    #plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmap = speed_cmap

    plotitem.pcolor_cmin = 0. 
    plotitem.pcolor_cmax = 20
    plotitem.add_colorbar = True
    plotitem.colorbar_label = 'knots'
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]
    plotitem.amr_patchedges_color = ['m','g','w']

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.show = False
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]


    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(-75,75,10)
    #plotitem.contour_nlevels = 10
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1,1,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-100,100]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'
    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        legend(('surface','topography'),loc='lower left')
        plot(t, 0*t, 'k')

    plotaxes.afteraxes = add_zeroline

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = range(40)
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once
    plotdata.html_movie_width = 700         # width used in JSAnimation

    return plotdata

    
