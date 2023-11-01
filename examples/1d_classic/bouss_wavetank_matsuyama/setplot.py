

import os, sys

try:
    from clawpack.geoclaw_1d import geoplot
except:
    print('Could not import from geoclaw_1d')


from clawpack.geoclaw_1d.nonuniform_grid_tools import make_mapc2p
import numpy

try:
    fname = '_output/fort.hmax'
    d = numpy.loadtxt(fname)
    etamax = numpy.where(d[:,1]>1e-6, d[:,2], numpy.nan)
    xmax = d[:,0]
    jmax = numpy.where(d[:,1]>0)[0].max()
    print("run-in = %8.2f,  run-up = %8.2f" % (d[jmax,0],d[jmax,2]))
    print('Loaded hmax from ',fname)
except:
    xmax = None
    print("Failed to load fort.hmax")

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


    def fixticks1(current_data):
        from pylab import ticklabel_format, grid,tight_layout
        ticklabel_format(useOffset=False)
        grid(True)
        tight_layout()
        #import pdb; pdb.set_trace()

    def fixticks(current_data):
        from pylab import ticklabel_format, plot,grid,gca
        ticklabel_format(useOffset=False)
        if xmax is not None:
            plot(xmax, etamax, 'r')
        grid(True)
        

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.3,0.6]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks

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
    plotaxes.title = 'Velocity'
    plotaxes.afteraxes = fixticks1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.velocity
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2 is not None:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.velocity
        plotitem.color = 'r'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-5,2]
    plotaxes.title = 'Full depth'
    plotaxes.afteraxes = fixticks1
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


    #----------

    plotfigure = plotdata.new_plotfigure(name='shore', figno=1)
    #plotfigure.kwargs = {'figsize':(9,11)}
    plotfigure.show = False
    

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [0,80e3]
    plotaxes.ylimits = [-4,4]
    plotaxes.title = 'Zoom on shelf'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.plot_var = geoplot.surface
    #plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    #plotaxes.xlimits = [-2000,2000]
    plotaxes.xlimits = [-1000,1000]
    #plotaxes.ylimits = [-10,40]
    plotaxes.ylimits = [-20,60]
    plotaxes.title = 'Zoom around shore'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
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
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.kwargs = {'figsize':(12,3)}
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = [-2,2]
    plotaxes.title = 'Surface elevation eta'

    def addgrid(current_data):
        from pylab import grid, xlim
        grid(True)
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
    plotaxes.afteraxes = addgrid

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2  # eta
    plotitem.plotstyle = 'b-'


    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output
    plotdata.parallel = True

    return plotdata

