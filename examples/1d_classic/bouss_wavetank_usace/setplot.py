
import os, sys

try:
    from clawpack.geoclaw_1d import geoplot
except:
    print('Could not import from geoclaw_1d')


from clawpack.geoclaw_1d.nonuniform_grid_tools import make_mapc2p
from clawpack.clawutil.data import ClawData
import numpy


if 0:
    try:
        fname = '_output/fgmax.txt'
        d = numpy.loadtxt(fname)
        etamax = numpy.where(d[:,1]>1e-6, d[:,3], numpy.nan)
        xmax = d[:,0]
        jmax = numpy.where(d[:,1]>0)[0].max()
        print("run-in = %8.2f,  run-up = %8.2f" % (d[jmax,0],d[jmax,3]))
        print('Loaded hmax from ',fname)
    except:
        xmax = None
        print("Failed to load fort.hmax")
else:
    xmax = None

xlimits = [-14,8.19]


outdir2 = None
#outdir2 = os.path.abspath('_output_ms')

if outdir2:
    if not os.path.isdir(outdir2):
        outdir2 = None

def setplot(plotdata):

    plotdata.clearfigures()

    outdir1 = plotdata.outdir
    mapc2p1, mx_edge, xp_edge = make_mapc2p(os.path.join(outdir1,'celledges.data'))
    
    from clawpack.amrclaw.data import GaugeData 
    setgauges = GaugeData() 
    setgauges.read(outdir1)
    gauge_xc = {}
    for k in range(len(setgauges.gauges)):
        gauge = setgauges.gauges[k]
        gaugeno = gauge[0]
        gauge_xc[gaugeno] = gauge[1]
    

    if outdir2:
        mapc2p2, mx_edge, xp_edge = make_mapc2p(os.path.join(outdir2,'celledges.data'))


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
        
    def velocity(current_data):
        from pylab import where
        q = current_data.q
        u = where(q[0,:]>1e-3, q[1,:] / q[0,:], 0.)
        return u

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.05,0.2]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.surface
        plotitem.plotstyle = 'r--'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p2

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.6,0.6]
    plotaxes.title = 'Velocity'
    plotaxes.afteraxes = fixticks1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.velocity
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.velocity
        plotitem.plotstyle = 'r--'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p2

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.25,0.1]
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
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    
    def fixgauge(current_data):
        from pylab import grid, title
        grid(True)
        gaugeno = current_data.gaugeno
        xc = gauge_xc[gaugeno]
        xp = mapc2p1(xc)
        print('+++ xc,xp:', xc,xp)
        title('Surface elevation at Gauge %i, x = %.0f m' \
              % (gaugeno, xp))

    plotaxes.afteraxes = fixgauge
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.02, 0.08]
    plotaxes.title = 'Surface elevation eta'
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

