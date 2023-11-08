
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import matplotlib.pyplot as plt
import os,sys

from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p


from clawpack.visclaw import gridtools  # would work with v5.7.0

# for comparing transects to 1d results
outdir_1d = os.path.abspath('1d_radial/_output')
#outdir_1d = None
print('Comparing to 1d solution in ', outdir_1d)

if outdir_1d is not None:
    fname_celledges = os.path.abspath('1d_radial/celledges.data')
    mapc2p1, mx_edge, xp_edge = make_mapc2p(fname_celledges)

xmin = 0.
xmax = 130e3



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
    plotitem.pcolor_cmin = -4.
    plotitem.pcolor_cmax = 4.
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1]
    plotitem.amr_patchedges_color = ['r','m','b','k']
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.show = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1]
    #plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffeeee']
    plotitem.amr_patch_bgcolor = [[1,.8,.8,.5],[.8,1,.8,.5]]
    plotitem.amr_patchedges_color = ['r','m','b','k']
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.contour_levels = [-4,-2,2,4]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]  
    

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo  #geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-2900,-110,-99,-1]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.schlieren_cmin = 0.0
    plotitem.schlieren_cmax = 2.0

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
        
        if 0:
            eps = 1.
            xout = linspace(xmin,xmax,1001)
            yout = zeros(xout.shape) + eps
            etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            hout = gridtools.grid_output_2d(framesoln, 0, xout, yout)
            zetaout = where(hout>0.001, etaout, nan)
            plot(xout, zetaout, 'b', label='along y=%g' % eps)
        
        if 0:
            eps = 20e3
            yout = linspace(xmin,xmax,1001)
            xout = zeros(yout.shape) + eps
            etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            hout = gridtools.grid_output_2d(framesoln, 0, xout, yout)
            zetaout = where(hout>0.001, etaout, nan)            
            plot(yout, pout, 'c', label='along x=%g' % eps)

        if 1:
            xout = linspace(xmin,xmax,14150)
            yout = xout
            etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            hout = gridtools.grid_output_2d(framesoln, 0, xout, yout)
            zetaout = where(hout>0.001, etaout, nan)
            plot(xout*sqrt(2), zetaout, 'm', label='along x=y')
            Bout = etaout - hout
            plot(xout*sqrt(2), Bout, 'g', label='topo along x=y')
            #yout = -xout
            #pout = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            #plot(abs(xout)*sqrt(2), pout, 'r', label='along x=-y')
            #plot(xout*sqrt(2), pout, 'r', label='along x=-y')

        if outdir_1d is not None:
            frameno = current_data.frameno
            fortt_file = '%s/fort.t%s' % (outdir_1d, str(frameno).zfill(4))
            with open(fortt_file,'r') as f:
                line = f.read()
                t1d = float(line.split()[0])
            if abs(t1d - current_data.t) < 1e-8:
                fortq_file = '%s/fort.q%s' % (outdir_1d, str(frameno).zfill(4))
                q1d = loadtxt(fortq_file, skiprows=5)
                mx1d = q1d.shape[0]
                dx1d = 1./mx1d
                xc1d = arange(dx1d/2, 1, dx1d)
                xp1d = mapc2p1(xc1d)
                h1d = q1d[:,0]
                eta1d = q1d[:,2]
                zeta1d = where(h1d > 0.001, eta1d, nan)
                plot(xp1d, zeta1d, 'k', label='1d radial')
            else:
                print('times do not match, not plotting 1d solution')
        legend()
        xlabel('radial distance')
        #xlim(0,xmax)
        xlim(xmin,xmax)
        ylim(-15,15)
        grid(True)
    plotaxes.afteraxes = plot_xsec



    #-----------------------------------------
    # Figure for surface and transect side by side
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='For paper', figno=120)
    plotfigure.figsize = (13,6)
    plotfigure.facecolor = 'w'
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    #plotaxes.scaled = True
    plotaxes.axescmd = 'axes([.1,.1,.35,.9])'

    axis_planview = { \
        18: (70e3,90e3,70e3,90e3),
        21: (77e3,86e3,77e3,86e3),
        23: (79e3,89e3,79e3,89e3),
        26: (80e3,90e3,80e3,90e3)
        }

    def fix_planview(current_data):
        from pylab import plot, title, xticks, axis
        frameno = current_data.frameno
        plot([0,124e3], [0,124e3], 'k--', linewidth=0.7)
        xticks(rotation=20)
        axis(axis_planview.get(frameno, 'scaled'))
        title('Eta at t = %.0f' % current_data.t, fontsize=15)
    plotaxes.afteraxes = fix_planview

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.show = False
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -4.
    plotitem.pcolor_cmax = 4.
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.5
    plotitem.colorbar_label = 'meters'
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,1,1,1]
    plotitem.amr_patchedges_color = ['r','m','b','yellow','k']
    #plotitem.amr_data_show = [0,0,1,1]  # for paper times > 1800
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo  #geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-2900,-110,-99,-1]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # cross section compared to 1d_radial
    #-----------------------------------------
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.48,.25,.45,.6])'
    plotaxes.title = '2D Transect / 1D Solution'
    #plotaxes.scaled = True

    axis_transect = { \
        6: (0,124e3,-5,5),
        12: (0,124e3,-5,5),
        18: (104e3,121e3,-4,6),
        21: (110e3,122e3,-4,9),
        23: (114e3,122e3,-4,8),
        26: (114e3,122e3,-4,8),
        }


    def plot_xsec(current_data):
        from pylab import plot,linspace,zeros,ones,legend,xlabel,sqrt,grid,xlim
        from pylab import nan,where,ylim,loadtxt,arange,ylabel,gca,axis,savefig
        from pylab import xticks,title
        from clawpack.pyclaw import Solution
        t = current_data.t
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)

        if 1:
            #xout = linspace(xmin,xmax,14150)
            xout = linspace(xmin,xmax,50000)
            yout = xout
            etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout)
            hout = gridtools.grid_output_2d(framesoln, 0, xout, yout)
            zetaout = where(hout>0.001, etaout, nan)
            plot(xout*sqrt(2), zetaout, 'm', label='2D along x=y')
            Bout = etaout - hout
            plot(xout*sqrt(2), Bout, 'g', label='topo along x=y')

        if outdir_1d is not None:
            frameno = current_data.frameno
            fortt_file = '%s/fort.t%s' % (outdir_1d, str(frameno).zfill(4))
            with open(fortt_file,'r') as f:
                line = f.read()
                t1d = float(line.split()[0])
            if abs(t1d - current_data.t) < 1e-8:
                fortq_file = '%s/fort.q%s' % (outdir_1d, str(frameno).zfill(4))
                q1d = loadtxt(fortq_file, skiprows=5)
                mx1d = q1d.shape[0]
                dx1d = 1./mx1d
                xc1d = arange(dx1d/2, 1, dx1d)
                xp1d = mapc2p1(xc1d)
                h1d = q1d[:,0]
                eta1d = q1d[:,2]
                zeta1d = where(h1d > 0.001, eta1d, nan)
                plot(xp1d, zeta1d, 'k', label='1D radial')
            else:
                print('times do not match, not plotting 1d solution')
                
        legend(loc='upper left', framealpha=1, fontsize=12)
        xlabel('radial distance', fontsize=12)
        xticks(rotation=20)
        ax = gca()
        ax.set_yticklabels([])
        ax2 = ax.secondary_yaxis('right')
        ax2.set_ylabel('meters', fontsize=12)
        #xlim(0,xmax)
        xlim(xmin,xmax)
        if 0:
            if t == 0.:
                ylim(-120,30)
            else:
                ylim(-15,15)
        frameno = current_data.frameno
        axis(axis_transect.get(frameno, (0,121e3,-15,15)))
        title('Eta on transect at t = %.0f' % current_data.t, fontsize=15)
        grid(True)
        if frameno in axis_transect.keys():
            fname = 'radocean%s.pdf' % str(frameno).zfill(2)
            savefig(fname, bbox_inches='tight')
            print('Created ',fname)

    plotaxes.afteraxes = plot_xsec



    #-----------------------------------------
    # Figure for radial velocity
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Radial momentum', figno=21)
    #plotfigure.show = False
    plotfigure.figsize = (8,8)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Radial momentum'
    plotaxes.scaled = True
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    def radial_momentum(current_data):
        from numpy import ma, sqrt
        q = current_data.q
        #h = q[0,:,:]
        hu = q[1,:,:]
        hv = q[2,:,:]
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        r = ma.masked_where(abs(r) < 1e-8, r)
        hs = (hu*x + hv*y) / r
        return hs
        
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = radial_momentum
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -800
    plotitem.pcolor_cmax = 800
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1


    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-990, -800, -600, -400, -210]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for cross section compared to 1d_radial
    #-----------------------------------------

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transects with 1d solution (black)'
    #plotaxes.scaled = True

    def radial_speed(q):
        hu = q[1,:,:]
        hv = q[2,:,:]
        hs = sqrt(hu**2 + hv**2)
        return hs
        
    def plot_xsec(current_data):
        from pylab import plot,linspace,zeros,ones,legend,xlabel,sqrt,grid,xlim
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3
        xout = linspace(xmin,xmax,1001)
        yout = zeros(xout.shape) + eps

        h = gridtools.grid_output_2d(framesoln, 0, xout, yout)
        eta = gridtools.grid_output_2d(framesoln, 5, xout, yout)
        B = eta-h
        plot(xout, B, 'g', linewidth=2, label='topo B on y=0')

        hu = gridtools.grid_output_2d(framesoln, 1, xout, yout)
        hv = gridtools.grid_output_2d(framesoln, 2, xout, yout)
        hs = (hu*xout + hv*yout) / sqrt(xout**2 + yout**2)
        plot(xout, hs, 'b', label='along y=0')
        yout = linspace(xmin,xmax,1001)
        xout = zeros(yout.shape) + eps
        hu = gridtools.grid_output_2d(framesoln, 1, xout, yout)
        hv = gridtools.grid_output_2d(framesoln, 2, xout, yout)
        hs = (hu*xout + hv*yout) / sqrt(xout**2 + yout**2)
        plot(yout, hs, 'c', label='along x=0')

        if 0:
            xout = linspace(xmin,xmax,1415)
            yout = xout
            hu = gridtools.grid_output_2d(framesoln, 1, xout, yout)
            hv = gridtools.grid_output_2d(framesoln, 2, xout, yout)
            hs = (hu*xout + hv*yout) / sqrt(xout**2 + yout**2)
            plot(xout*sqrt(2), hs, 'm', label='along x=y')
            yout = -xout
            hu = gridtools.grid_output_2d(framesoln, 1, xout, yout)
            hv = gridtools.grid_output_2d(framesoln, 2, xout, yout)
            hs = (hu*xout + hv*yout) / sqrt(xout**2 + yout**2)
            plot(xout*sqrt(2), hs, 'r', label='along x=-y')

        legend()
        xlabel('radial distance')
        xlim(xmin,xmax)
        grid(True)
    plotaxes.afteraxes = plot_xsec



    #-----------------------------------------
    # Figure for huc corrections
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='huc', figno=25)
    #plotfigure.show = False
    plotfigure.figsize = (8,8)
    plotfigure.facecolor = 'w'
    ihuc = 3

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'hu correction'
    plotaxes.scaled = True
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = ihuc
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.002
    plotitem.pcolor_cmax =  0.003  
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-990, -800, -600, -400, -210]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for cross section 
    #-----------------------------------------
    #plotfigure = plotdata.new_plotfigure(name='Radial section', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transects'
    #plotaxes.scaled = True

    def plot_xsec(current_data):
        from pylab import plot,linspace,zeros,ones,legend,xlabel,sqrt,grid,xlim
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3
        xout = linspace(xmin,xmax,1001)
        yout = zeros(xout.shape) + eps
        pout = gridtools.grid_output_2d(framesoln, ihuc, xout, yout)
        #plot(abs(xout), pout, 'b', label='along y=0')
        plot(xout, pout, 'b', label='along y=0')
        yout = linspace(xmin,xmax,1001)
        xout = zeros(yout.shape) + eps
        pout = gridtools.grid_output_2d(framesoln, ihuc, xout, yout)
        #plot(abs(yout), pout, 'c', label='along x=0')
        plot(yout, pout, 'c', label='along x=0')
        xout = linspace(xmin,xmax,1415)
        yout = xout
        pout = gridtools.grid_output_2d(framesoln, ihuc, xout, yout)
        #plot(abs(xout)*sqrt(2), pout, 'm', label='along x=y')
        #plot(xout*sqrt(2), pout, 'm', label='along x=y')
        yout = -xout
        pout = gridtools.grid_output_2d(framesoln, ihuc, xout, yout)
        #plot(abs(xout)*sqrt(2), pout, 'r', label='along x=-y')
        #plot(xout*sqrt(2), pout, 'r', label='along x=-y')
        legend()
        xlabel('radial distance')
        #xlim(0,xmax)
        xlim(xmin,xmax)
        grid(True)
    plotaxes.afteraxes = plot_xsec


    #-----------------------------------------
    # Figure for hvc corrections
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='hvc', figno=26)
    plotfigure.show = False
    plotfigure.figsize = (8,8)
    plotfigure.facecolor = 'w'
    ihvc = 4

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'hv correction'
    plotaxes.scaled = True
    plotaxes.axescmd = 'axes([.15,.5,.7,.45])'

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = ihvc
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.002
    plotitem.pcolor_cmax =  0.003  
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [-990, -800, -600, -400, -210]
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]  
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for cross section 
    #-----------------------------------------
    #plotfigure = plotdata.new_plotfigure(name='Radial section', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('radial slice')
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'Transects'
    #plotaxes.scaled = True

    def plot_xsec(current_data):
        from pylab import plot,linspace,zeros,ones,legend,xlabel,sqrt,grid,xlim
        from clawpack.pyclaw import Solution
        pd = current_data.plotdata
        frameno = current_data.frameno
        framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
        eps = 1e-3
        xout = linspace(xmin,xmax,1001)
        yout = zeros(xout.shape) + eps
        pout = gridtools.grid_output_2d(framesoln, ihvc, xout, yout)
        #plot(abs(xout), pout, 'b', label='along y=0')
        plot(xout, pout, 'b', label='along y=0')
        yout = linspace(xmin,xmax,1001)
        xout = zeros(yout.shape) + eps
        pout = gridtools.grid_output_2d(framesoln, ihvc, xout, yout)
        #plot(abs(yout), pout, 'c', label='along x=0')
        plot(yout, pout, 'c', label='along x=0')
        xout = linspace(xmin,xmax,1415)
        yout = xout
        pout = gridtools.grid_output_2d(framesoln, ihvc, xout, yout)
        #plot(abs(xout)*sqrt(2), pout, 'm', label='along x=y')
        #plot(xout*sqrt(2), pout, 'm', label='along x=y')
        yout = -xout
        pout = gridtools.grid_output_2d(framesoln, ihvc, xout, yout)
        #plot(abs(xout)*sqrt(2), pout, 'r', label='along x=-y')
        #plot(xout*sqrt(2), pout, 'r', label='along x=-y')
        legend()
        xlabel('radial distance')
        #xlim(0,xmax)
        xlim(xmin,xmax)
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
