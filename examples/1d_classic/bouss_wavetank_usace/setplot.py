
import os, sys
from clawpack.visclaw import geoplot
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
from clawpack.clawutil.data import ClawData
import numpy


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


    def fix_layout(current_data):
        from pylab import tight_layout
        tight_layout()

        
    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.05,0.2]
    plotaxes.grid = True
    plotaxes.title = 'Surface displacement'

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
    plotaxes.grid = True
    plotaxes.title = 'Velocity'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.u_velocity
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p1

    if outdir2:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.u_velocity
        plotitem.plotstyle = 'r--'
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p2

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-0.25,0.1]
    plotaxes.grid = True
    plotaxes.title = 'Full depth'
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
    plotfigure.clf_each_gauge = True
    plotfigure.figsize = (12,3)

    plotaxes = plotfigure.new_plotaxes()
    
    def fixgauge(current_data):
        from pylab import title
        gaugeno = current_data.gaugeno
        xc = gauge_xc[gaugeno]
        xp = mapc2p1(xc)
        #print('+++ xc,xp:', xc,xp)
        title('Surface elevation at Gauge %i, x = %.0f m' \
              % (gaugeno, xp))

    plotaxes.afteraxes = fixgauge
    plotaxes.xlimits = [0,25]
    plotaxes.ylimits = [-0.02, 0.08]
    plotaxes.grid = True
    plotaxes.title = 'Surface elevation eta'
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
