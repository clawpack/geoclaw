

import os, sys

from clawpack.visclaw import geoplot
import numpy
from scipy import fft
from numpy import exp,cos,pi
from pylab import figure,clf,plot,title

xlimits = [0, 5e3]

#outdir2 = '_output_100m_sgn'  # second outdir for comparison plots, if desired
outdir2 = None


def setplot(plotdata):

    plotdata.clearfigures()

    def set_dry_tol(current_data):
        current_data.user['dry_tolerance'] = 1e-6

    plotdata.beforeframe = set_dry_tol

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    #plotaxes.ylimits = [-5,15]
    plotaxes.title = 'Surface displacement'
    plotaxes.grid = True
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'

    if outdir2 is not None:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.surface
        plotitem.color = 'r'
        #plotitem.plotstyle = '-+'


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    #plotaxes.ylimits = [-1,1]
    plotaxes.title = 'Velocity'
    plotaxes.grid = True
        

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.u_velocity
    plotitem.color = 'b'

    if outdir2 is not None:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.outdir = outdir2
        plotitem.plot_var = geoplot.u_velocity
        plotitem.color = 'r'



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
