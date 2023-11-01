

import os, sys
from clawpack.visclaw import geoplot
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
import numpy


fname_celledges = os.path.abspath('celledges.data')


xlimits = [-300,150]

def setplot(plotdata=None):

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()

    mapc2p, mx_edge, xp_edge = make_mapc2p(fname_celledges)

    def mapc2p_km(xc):
        x_m = mapc2p(xc)
        x_km = x_m / 1000.   # convert to km
        return x_km


    def add_annotations(current_data):
        from pylab import ticklabel_format, plot,grid,ones,sqrt, \
            legend,title,ylabel,text
        ticklabel_format(useOffset=False)

        hl = 3200.
        hr = 200.
        greens = (hl/hr)**(0.25)
        print('greens = ',greens)
        #plot(current_data.x, greens*ones(current_data.x.shape),'g--')
        plot(xlimits,[greens,greens],'g--', label='$C_g$, Greens Law')
        ctrans = 2*sqrt(hl)/(sqrt(hl)+sqrt(hr))
        crefl = (sqrt(hl)-sqrt(hr))/(sqrt(hl)+sqrt(hr))
        print('ctrans = ',ctrans)
        plot(xlimits,[ctrans,ctrans],'r--', label='$C_T$, Transmission coefficient')
        print('crefl = ',crefl)
        plot(xlimits,[crefl,crefl],'m--', label='$C_R$, Reflection coefficient')
        legend(loc='upper left')
        title('')
        ylabel('meters', fontsize=14)
        if current_data.frameno == 0:
            text(-95,-0.4,'$\longrightarrow$',fontsize=20)
            text(-95,-0.6,'Incident')
        h = current_data.q[0,:]
        mx2 = int(round(len(h)/2.))
        etamax2 = (h[:mx2] - hl).max()
        print('mx2 = %i, etamax2 = %g' % (mx2,etamax2))
        if (current_data.frameno == 5) and (etamax2 > 0.1):
            text(-190,-0.5,'$\longleftarrow$',fontsize=20)
            text(-190,-0.7,'Reflected')
            text(30,-0.5,'$\longrightarrow$',fontsize=20)
            text(15,-0.7,'Transmitted')
        if (current_data.frameno == 6) and (etamax2 > 0.1):
            text(-260,-0.5,'$\longleftarrow$',fontsize=20)
            text(-260,-0.7,'Reflected')
            text(40,-0.5,'$\longrightarrow$',fontsize=20)
            text(25,-0.7,'Transmitted')
        elif (current_data.frameno == 6):
            text(-20,-0.5,'$\longleftarrow$',fontsize=20)
            text(-20,-0.7,'Reflected')
            text(70,-0.5,'$\longrightarrow$',fontsize=20)
            text(65,-0.7,'Transmitted')
        

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize':(7,6.5)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.4,.8,.5])' #'subplot(211)'
    plotaxes.xlimits = xlimits
    #plotaxes.xlimits = [-100e3,-20e3]
    plotaxes.ylimits = [-1,3]
    plotaxes.title = 'Surface displacement'
    plotaxes.afteraxes = add_annotations

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km


    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = 1
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.1,.8,.2])' #'subplot(212)'
    plotaxes.xlimits = xlimits

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'g'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km

    def fix_topo_plot(current_data):
        from pylab import title,xlabel
        title('')
        xlabel('kilometers', fontsize=14)
    plotaxes.afteraxes = fix_topo_plot

    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p_km



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                                         type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Eta'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2
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
