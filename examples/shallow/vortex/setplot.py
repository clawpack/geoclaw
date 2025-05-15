#!/usr/bin/env python

import os

import numpy as np
import matplotlib.pyplot as plt

# import clawpack.visclaw.colormaps as colormap
import clawpack.visclaw.gaugetools as gaugetools
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata
import clawpack.visclaw.geoplot as geoplot
import clawpack.geoclaw.surge.plot as surgeplot

try:
    from setplotfg import setplotfg
except:
    setplotfg = None

M = .5
g = 1.0
c1 = 0.04
c2 = 0.02
alpha = np.pi / 6.
x0 = -20.
y0 = -10.
def exact_solution(x, y, t):
    f = lambda x,y,t: -c2*((x-x0-M*t*np.cos(alpha))**2+(y-y0-M*t*np.sin(alpha))**2)
    h = lambda x,y,t: 1.-c1**2/(4.*c2*g)*np.exp(2.*f(x,y,t))
    u = lambda x,y,t: M*np.cos(alpha)+c1*(y-y0-M*t*np.sin(alpha))*np.exp(f(x,y,t))
    v = lambda x,y,t: M*np.sin(alpha)-c1*(x-x0-M*t*np.cos(alpha))*np.exp(f(x,y,t))

    return h(x, y, t), u(x, y, t), v(x, y, t)

def exact_vorticity(x, y, t):
    f = lambda x,y,t: -c2*((x-x0-M*t*np.cos(alpha))**2+(y-y0-M*t*np.sin(alpha))**2)
    f_x = lambda x, y, t: -2 * c2 * (x - x0 - M * t * np.cos(alpha))
    f_y = lambda x, y, t: -2 * c2 * (y - y0 - M * t * np.sin(alpha))

    u_y = lambda x, y, t:  c1 * np.exp(f(x, y, t)) * (1 + (y - y0 - M * t * np.sin(alpha)) * f_y(x, y, t))
    v_x = lambda x, y, t: -c1 * np.exp(f(x, y, t)) * (1 + (x - x0 - M * t * np.cos(alpha)) * f_x(x, y, t))

    return v_x(x, y, t) - u_y(x, y, t)

def extract_eta(cd):
    return cd.q[0, :, :] - 1

def eta_error(cd):
    h, u, v = exact_solution(cd.x, cd.y, cd.t)
    return extract_eta(cd) - (h - 1.0)

def speed_error(cd):
    speed = surgeplot.water_speed(cd)
    h, u, v = exact_solution(cd.x, cd.y, cd.t)
    return speed - np.sqrt(u**2 + v**2)

def vorticity_error(cd):
    omega = surgeplot.water_vorticity(cd)
    # h, u, v = exact_solution(cd.x, cd.y, cd.t)
    # delta = [cd.x[1, 0] - cd.x[0, 0], cd.y[0, 1] - cd.y[0, 0]]
    # exact_omega = surgeplot.extract_vorticity(delta, u, v)
    exact_omega = exact_vorticity(cd.x, cd.y, cd.t)
    return omega - exact_omega

# Setplot
def setplot(plotdata=None):
    """"""

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    # clear any old figures,axes,items data
    plotdata.clearfigures()
    plotdata.format = 'ascii'

    # Load data from output
    try:
        clawdata = clawutil.ClawInputData(2)
        clawdata.read(os.path.join(plotdata.outdir, 'claw.data'))
        physics = geodata.GeoClawData()
        physics.read(os.path.join(plotdata.outdir, 'geoclaw.data'))
        surge_data = geodata.SurgeData()
        surge_data.read(os.path.join(plotdata.outdir, 'surge.data'))
        friction_data = geodata.FrictionData()
        friction_data.read(os.path.join(plotdata.outdir, 'friction.data'))

        domain_limits = ((clawdata.lower[0], clawdata.upper[0]), 
                         (clawdata.lower[1], clawdata.upper[1]))
    except:
        # Assume that we are running PyClaw
        domain_limits = (None, None)    

    # Color limits
    surface_limits = [-0.02, 0.02]
    speed_limits = [0.3, 0.6]
    vorticity_limits = [-0.08, 0.08]
    eta_error_limits = [-0.02, 0.02]
    speed_error_limits = [-0.07, 0.07]
    vorticity_error_limits = [-0.007, 0.007]

    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    # Surface Figure
    plotfigure = plotdata.new_plotfigure(name="Surface")
    plotfigure.kwargs = {"figsize": (6.4, 4.8)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Surface"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]

    surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
    plotaxes.plotitem_dict['surface'].plot_var = extract_eta
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10

    # Speed Figure
    plotfigure = plotdata.new_plotfigure(name="Currents")
    plotfigure.kwargs = {"figsize": (6.4, 4.8)}
    plotfigure.show = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Currents"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]

    surgeplot.add_speed(plotaxes, bounds=speed_limits)
    plotaxes.plotitem_dict['speed'].amr_patchedges_show = [0] * 10

    # Vorticity
    plotfigure = plotdata.new_plotfigure(name="Vorticity")
    plotfigure.kwargs = {}
    plotfigure.show = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Vorticity"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]

    surgeplot.add_vorticity(plotaxes, bounds=vorticity_limits)
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10

    # ========================================================================
    # Error Plots
    def add_eta_error_title(cd, order=1):
        delta = (cd.x[1, 0] - cd.x[0, 0]) * (cd.y[0, 1] - cd.y[0, 0])
        error = np.linalg.norm(eta_error(cd), ord=order) * delta
        plt.gca().set_title(r"$\eta$ Error (t = {}), ".format(cd.t) + 
                            r"$|| E ||_{\ell_1} = $" +
                            f"{error.round(6)}")

    def add_speed_error_title(cd, order=1):
        delta = (cd.x[1, 0] - cd.x[0, 0]) * (cd.y[0, 1] - cd.y[0, 0])
        error = np.linalg.norm(speed_error(cd), ord=order) * delta
        plt.gca().set_title(r"Speed Error (t = {}), ".format(cd.t) + 
                            r"$|| E ||_{\ell_1} = $" +
                            f"{error.round(6)}")

    def add_vorticity_error_title(cd, order=1):
        delta = (cd.x[1, 0] - cd.x[0, 0]) * (cd.y[0, 1] - cd.y[0, 0])
        error = np.linalg.norm(vorticity_error(cd), ord=order) * delta
        plt.gca().set_title(r"$\omega$ Error (t = {}), ".format(cd.t) + 
                            r"$|| E ||_{\ell_1} = $" +
                            f"{error.round(6)}")

    # Surface
    plotfigure = plotdata.new_plotfigure(name="Surface Error")
    plotfigure.kwargs = {"figsize": (6.4, 4.8)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = r"$\eta$ Error"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]
    plotaxes.afteraxes = add_eta_error_title
    
    plotitem = plotaxes.new_plotitem(name='surface error', plot_type='2d_pcolor')
    plotitem.plot_var = eta_error
    plotitem.pcolor_cmax = eta_error_limits[1]
    plotitem.pcolor_cmin = eta_error_limits[0]
    plotitem.pcolor_cmap = plt.get_cmap("RdBu")
    plotitem.add_colorbar = True
    plotitem.colorbar_label = r"$\eta$ Error"

    # Speed
    plotfigure = plotdata.new_plotfigure(name="Speed Error")
    plotfigure.kwargs = {"figsize": (6.4, 4.8)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = r"Speed Error"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]
    plotaxes.afteraxes = add_speed_error_title
    
    plotitem = plotaxes.new_plotitem(name='speed error', plot_type='2d_pcolor')
    plotitem.plot_var = speed_error
    plotitem.pcolor_cmax = speed_error_limits[1]
    plotitem.pcolor_cmin = speed_error_limits[0]
    plotitem.pcolor_cmap = plt.get_cmap("RdBu")
    plotitem.add_colorbar = True
    plotitem.colorbar_label = r"Speed Error"

    # Vorticity
    plotfigure = plotdata.new_plotfigure(name="Vorticity Error")
    plotfigure.kwargs = {"figsize": (6.4, 4.8)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = r"Vorticity Error"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]
    plotaxes.afteraxes = add_vorticity_error_title
    
    plotitem = plotaxes.new_plotitem(name='vorticity error', plot_type='2d_pcolor')
    plotitem.plot_var = vorticity_error
    plotitem.pcolor_cmax = vorticity_error_limits[1]
    plotitem.pcolor_cmin = vorticity_error_limits[0]
    plotitem.pcolor_cmap = plt.get_cmap("RdBu")
    plotitem.add_colorbar = True
    plotitem.colorbar_label = r"Vorticity Error"

    # ========================================================================
    # Transects
    def compute_max(current_data, field=3, title=r"Field {} - $\max = {}$"):
        max_value = np.max(np.abs(current_data.q[field, :, :]))
        # plt.gca().set_title(title.format(field, max_value))

    def transect_eta(cd, y0=0.0):
        y = cd.y
        dy = cd.dy
        index = np.where(abs(y - y0) <= dy / 2.0)[1][0]
        if cd.q.shape[0] > 3:
            eta = surgeplot.extract_eta(cd.q[0, :, index], cd.q[3, :, index])
        else:
            eta = cd.q[0, :, index] - 1
        return cd.x[:, index], eta

    def transect_velocity(current_data, y0=0.0):
        y = current_data.y
        dy = current_data.dy
        index = np.where(abs(y - y0) <= dy / 2.0)[1][0]
        h = current_data.q[0, :, index]
        hu = current_data.q[1, :, index]
        hv = current_data.q[2, :, index]
        u = np.where(h > 1e-3, hu / h, np.zeros(h.shape))
        v = np.where(h > 1e-3, hv / h, np.zeros(h.shape))
        return current_data.x[:, index], u
        # return current_data.x[:, index], v

    # === Depth ===
    plotfigure = plotdata.new_plotfigure(name="Depth Transect")
    plotfigure.show = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Depth Transect"
    plotaxes.xlabel = "x (m)"
    plotaxes.ylabel = r"$h$"
    plotaxes.xlimits = domain_limits[0]
    # plotaxes.ylimits = [0.97, 1.01]
    plotaxes.ylimits = surface_limits
    plotaxes.grid = True
    # plotaxes.afteraxes = lambda cd: compute_max(cd, field=0)

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = transect_eta
    plotitem.plotstyle = 'ko-'
    plotitem.kwargs = {"markersize": 3}

    # === Momentum/Velocity ===
    plotfigure = plotdata.new_plotfigure(name="Velocity Transect")
    plotfigure.show = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Velocity Transect"
    plotaxes.xlabel = "x (m)"
    plotaxes.ylabel = r"$u$"
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = speed_limits
    plotaxes.grid = True
    # plotaxes.afteraxes = lambda cd: compute_max(cd, field=1)


    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = transect_velocity
    plotitem.plotstyle = 'bx-'
    plotitem.kwargs = {"markersize": 3}

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Gauge Surfaces', figno=300,
                                         type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.grid = True
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = surface_limits
    plotaxes.title = "Surface"
    plotaxes.ylabel = "Surface (m)"
    plotaxes.time_label = "t (s)"

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = surgeplot.gauge_surface

    #  Gauge Location Plot
    def gauge_location_afteraxes(cd):
        gaugetools.plot_gauge_locations(cd.plotdata, gaugenos='all',
                                        format_string='kx', fontsize=10,
                                        add_labels=True)

    plotfigure = plotdata.new_plotfigure(name="Gauge Locations")
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Gauge Locations'
    plotaxes.scaled = True
    plotaxes.xlimits = domain_limits[0]
    plotaxes.ylimits = domain_limits[1]
    plotaxes.afteraxes = gauge_location_afteraxes
    surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
    plotaxes.plotitem_dict['surface'].plot_var = extract_eta
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10

    # -----------------------------------------
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'        # list of gauges to print
    # plotdata.print_gaugenos = [10, 21, 32]   # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # parallel plotting

    return plotdata
