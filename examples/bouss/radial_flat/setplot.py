
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import matplotlib.pyplot as plt
import os,sys


from clawpack.visclaw import gridtools  # would work with v5.7.0

# for comparing transects to 1d results
outdir_1d = os.path.abspath('1d_radial/_output')
#outdir_1d = None
print('Comparing to 1d solution in ', outdir_1d)


xmin = 0.
xmax = 5000.



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
    #plotdata.format = 'binary'    # 'ascii' or 'binary' to match setrun.py


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=20)
    #plotfigure.show = False
    plotfigure.figsize = (8,8)
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1]
    plotitem.amr_patchedges_color = ['r','yellow','w','k']
    

    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------
    #plotfigure = plotdata.new_plotfigure(name='Radial section', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transects with 1d solution (black)'
    #plotaxes.scaled = True

    def plot_xsec(current_data):
        from pylab import plot,linspace,zeros,ones,legend,xlabel,sqrt,grid,xlim
        from pylab import nan,where,ylim,loadtxt,arange
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        

        if 1:
            xout = linspace(xmin,xmax,14150)
            yout = xout
            etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            hout = gridtools.grid_output_2d(framesoln, 0, xout, yout)
            zetaout = where(hout>0.001, etaout, nan)
            plot(xout*sqrt(2), zetaout, 'm', label='along x=y')

        if outdir_1d is not None:
            frameno = current_data.frameno
            framesoln = Solution(frameno,path=outdir_1d, file_format='ascii')
            t1d = framesoln.t
            if abs(t1d - current_data.t) < 1e-8:
                state = framesoln.states[0]
                x1d = state.grid.c_centers[0]
                q1d = state.q
                h1d = q1d[0,:]
                eta1d = q1d[2,:]
                zeta1d = where(h1d > 0.001, eta1d, nan)
                plot(x1d, zeta1d, 'k', label='1d radial')
            else:
                print('times do not match, not plotting 1d solution')

        legend()
        xlabel('radial distance')
        #xlim(0,xmax)
        xlim(xmin,xmax)
        ylim(-0.5,0.5)
        grid(True)
    plotaxes.afteraxes = plot_xsec



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
