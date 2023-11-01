

import os, sys
from clawpack.visclaw import geoplot
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
import numpy


xlimits = [-160,40]

fname_celledges = os.path.abspath('celledges.data')

outdir2 = None  # set to compare different models

def setplot(plotdata):

    plotdata.clearfigures()

    def set_dry_tol(current_data):
        current_data.user['dry_tolerance'] = 1e-6

    plotdata.beforeframe = set_dry_tol

    outdir1 = plotdata.outdir
    mapc2p1, mx_edge, xp_edge = make_mapc2p(fname_celledges)


    def fix_layout(current_data):
        from pylab import tight_layout
        tight_layout()


    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.figsize = (8,7)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.3,0.6]
    plotaxes.grid = True
    plotaxes.title = 'Surface displacement'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    #plotitem.plotstyle = '-+'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2 is not None:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.surface
        plotitem.color = 'r'
        #plotitem.plotstyle = '-+'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-1,1]
    plotaxes.grid = True
    plotaxes.title = 'Velocity'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.u_velocity
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2 is not None:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.u_velocity
        plotitem.color = 'r'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-5,2]
    plotaxes.grid = True
    plotaxes.title = 'Full depth'
    plotaxes.grid = True
    plotaxes.afteraxes = fix_layout
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [.5,.5,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.figsize = (12,3)
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = [-2,2]
    plotaxes.grid = True
    plotaxes.title = 'Surface elevation eta'

    def setxlim(current_data):
        from pylab import xlim
        gno = current_data.gaugeno
        if gno == 150:
            xlim(0,60)
        elif gno == 80:
            xlim(10,60)
        elif gno == 60:
            xlim(20,70)
        elif gno == 40:
            xlim(40,90)
        elif gno == 30:
            xlim(40,90)
        elif gno == 20:
            xlim(50,100)
        elif gno == 10:
            xlim(60,110)
    plotaxes.afteraxes = setxlim

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2  # eta
    plotitem.plotstyle = 'b-'


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
