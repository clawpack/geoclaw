

import os, sys
import numpy
from clawpack.visclaw import geoplot
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p

fname_celledges = os.path.abspath('celledges.data')

fname = '_output/fgmax.txt'
try:
    d = numpy.loadtxt(fname)
    xmax = d[:,0]
    Bmax = d[:,1]
    hmax = d[:,2]
    etamax = numpy.where(hmax>1e-3, hmax+Bmax, numpy.nan)
    jmax = numpy.where(hmax>1e-3)[0].max()
    print("run-in = %8.2f m,  run-up = %8.2f m" % (xmax[jmax],etamax[jmax]))
    print('Loaded hmax from ',fname)
except:
    xmax = None
    print('Failed to load ',fname)

#xmax = None # to suppress plotting max elevation as red curve 

xlimits = [-200e3,1e3]

def setplot(plotdata):

    plotdata.clearfigures()

    fname1 = os.path.join(plotdata.outdir, fname_celledges)
    mapc2p1, mx_edge, xp_edge = make_mapc2p(fname1)


    def fix_layout(current_data):
        from pylab import tight_layout
        tight_layout()

    def add_etamax(current_data):
        from pylab import plot
        from clawpack.visclaw.legend_tools import add_legend
        if xmax is not None:
            plot(xmax, etamax, 'r')
            add_legend(['max eta over simulation','surface elevation eta'],
                   ['r','b'], framealpha=1)

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.figsize = (8,6)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-5,5]
    plotaxes.title = 'Surface displacement at time h:m:s'
    plotaxes.grid = True

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.title = 'Full depth at time h:m:s'
    plotaxes.grid = True
    plotaxes.afteraxes = fix_layout

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = [0.4,0.4,1] # lighter blue
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1


    #----------

    plotfigure = plotdata.new_plotfigure(name='shore', figno=1)
    plotfigure.figsize = (8,6)
    #plotfigure.show = False
    

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [-50e3,1e3]
    plotaxes.ylimits = [-20,20]
    plotaxes.title = 'Zoom on shelf at time h:m:s'
    plotaxes.grid = True
    plotaxes.afteraxes = add_etamax

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [-100,300]
    plotaxes.ylimits = [-10,15]
    plotaxes.title = 'Zoom around shore at time h:m:s'
    plotaxes.grid = True
    plotaxes.afteraxes = fix_layout

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='depth', figno=300, \
                    type='each_gauge')

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = 1/3600.  # convert seconds to hours
    plotaxes.time_label = 'time (hours) post-quake'
    plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = [-2.0, 2.0]
    plotaxes.title = 'Water depth'
    plotaxes.grid = True

    # Plot depth as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
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
