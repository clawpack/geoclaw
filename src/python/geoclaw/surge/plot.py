# encoding: utf-8
r"""
Plotting routines for storm surge simulations with GeoClaw

:Authors:
    Kyle Mandli (2012-10-23) Initial version
"""
# ==============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD)
#  license
#                     http://www.opensource.org/licenses/
# ==============================================================================

from __future__ import absolute_import
from __future__ import print_function

import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as mlines
import pandas

import clawpack.visclaw.colormaps as colormaps
import clawpack.visclaw.gaugetools as gaugetools
import clawpack.visclaw.geoplot as geoplot
import clawpack.geoclaw.data as geodata
# import clawpack.geoclaw.surge.storm

# TODO:  Assign these based on data files
bathy_index = 0
friction_field = 3
wind_field = 4
pressure_field = 6

surface_cmap = plt.get_cmap("bwr")
speed_cmap = plt.get_cmap('PuBu')
friction_cmap = plt.get_cmap('YlOrRd')
velocity_cmap = plt.get_cmap('PiYG')
vorticity_cmap = plt.get_cmap('PRGn')
wind_cmap = plt.get_cmap('PuBu')
pressure_cmap = plt.get_cmap('PuBu')
land_cmap = geoplot.land_colors

surge_data = geodata.SurgeData()

# ===========================
#  General Drawing Functions
# ===========================
def draw_box(ax, box, style='r-', fill=True):
    r"""Draw the box specified

    :Input:
     - *ax* (matplotlib.Axes) axes object to plot on
     - *box* (list) set of coordinates defined box.  First two coordinates are 
       the lower left corner and the last two coordiantes are the upper corner.
     - *style* (str) string that will be used to specify the style of the lines.
       Default = *'r-'*.
     - *fill* (bool) If True draws a filled box.  Not implemented.

    :TODO:
     - Handle fill request
    """
    ax.plot([box[0], box[2]], [box[1], box[1]], style) # Bottom
    ax.plot([box[0], box[2]], [box[3], box[3]], style) # Top
    ax.plot([box[0], box[0]], [box[1], box[3]], style) # Left
    ax.plot([box[2], box[2]], [box[1], box[3]], style) # Right


# ==============================
#  Track Plotting Functionality
# ==============================
class track_data(object):
    """Read in storm track data from run output"""

    def __init__(self, path=None):
        if path is None:
            path = "fort.track"

        try:
            self._path = path
            self._data = np.loadtxt(self._path)
        except:
            self._data = None

    def get_track(self, frame):
        """Return storm location for frame requested"""

        # If data was not load successfully return None
        if self._data is None or len(self._data.shape) < 2:
            return None, None, None

        # If it appears that our data is not long enough, try reloading file
        if self._data.shape[0] < frame + 1:
            self._data = np.loadtxt(self._path)

            # Check to make sure that this fixed the problem
            if self._data.shape[0] < frame + 1:
                warnings.warn(f" *** WARNING *** Could not find track data",
                               " for frame {frame}.")
                return None, None, None

        return self._data[frame, 1:]


# ==========================================================================
# Gauge functions
# ==========================================================================
def gauge_locations(current_data, gaugenos='all'):
    gaugetools.plot_gauge_locations(current_data.plotdata,
                                    gaugenos=gaugenos, format_string='kx',
                                    add_labels=True, xoffset=0.02,
                                    yoffset=0.02)


def gauge_dry_regions(cd, dry_tolerance=1e-16):
    """Masked array of zeros where gauge is dry."""
    return np.ma.masked_where(np.abs(cd.gaugesoln.q[0, :]) > dry_tolerance,
                              np.zeros(cd.gaugesoln.q[0, :].shape))


def gauge_surface(cd, dry_tolerance=1e-16):
    """Sea surface at gauge masked when dry."""
    return np.ma.masked_where(np.abs(cd.gaugesoln.q[0, :]) < dry_tolerance,
                              cd.gaugesoln.q[3, :])


def plot_landfall_gauge(gauge, axes, landfall=0.0, style='b', kwargs={}):
    """Plot gauge data on the axes provided

    This will transform the plot so that it is relative to the landfall value
    provided.

    This can be done using `plotaxes.time_scale` instead so this function will
    be deprecated and removed in a future release.
    """
    axes = plt.gca()

    # Add GeoClaw gauge data
    t = sec2days(gauge.t - landfall)
    axes.plot(t, gauge.q[3, :], style, **kwargs)


# ========================================================================
#  Surge related helper functions
# ========================================================================
def days_figure_title(cd, land_fall=0.0, new_time=False):
    r"""Helper function that puts the time relative to landfall in title

    New version of title is available if *new_time = True*
    """
    if new_time:
        if cd.t < land_fall:
            sign = "-"
        else:
            sign = " "
        minutes, seconds = divmod(abs(np.round(cd.t - land_fall)), 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        days = int(days)
        hours = int(hours)
        minutes = int(minutes)
        if cd.t < 0:
            sign = "-"
        else:
            sign = " "
        title = cd.plotaxes.title
        plt.title(f'{title} at t = {sign}{days:d}, {hours:02d}:{minutes:02d}')
    else:
        t = (cd.t - land_fall) / (60**2 * 24)
        days = int(t)
        hours = (t - int(t)) * 24.0

        title = cd.plotaxes.title
        plt.title('%s at day %3i, hour %2.1f' % (title, days, hours))


def surge_afteraxes(current_data, track, land_fall=0.0, plot_direction=False,
                    plot_track=False, style='ro', track_style='k', 
                    new_time=False, kwargs={}):
    """Default surge plotting after axes function

    Includes changing the title to something relative to landfall and plotting
    the location of the storm eye according to the track object.  Optional 
    plotting includes the direction of the storm with *plot_direction* and 
    *plot_track* which plots the entire track of the storm.
    """
    track_data = track.get_track(current_data.frameno)

    if track_data[0] is not None and track_data[1] is not None:
        ax = plt.gca()
        ax.plot(track_data[0], track_data[1], style, **kwargs)
        # Plot direction
        if plot_direction:
            ax.quiver(track_data[0], track_data[1],
                        np.cos(track_data[2]), np.sin(track_data[2]))
        # Plot full track
        if plot_track:
            # TO DO: Could add categorization but really should use storm object
            # for this
            for i in range(1, track._data.shape[0]):
                ax.plot(track._data[i-1:i+1, 1], track._data[i-1:i+1, 2], 
                        track_style)
        
    days_figure_title(current_data, land_fall=land_fall, new_time=new_time)


def friction(cd):
    return cd.aux[friction_field, :, :]


def wind_x(cd):
    # print(cd.aux[wind_field, :, :])
    return cd.aux[wind_field, :, :]


def wind_y(cd):
    # print(cd.aux[wind_field+1, :, :])
    return cd.aux[wind_field+1, :, :]


def wind_speed(cd):
    return np.sqrt(wind_x(cd)**2 + wind_y(cd)**2)


def pressure(cd):
    # The division by 100.0 is to convert from Pa to millibars
    return cd.aux[pressure_field, :, :] / 100.0


def storm_radius(cd, track):
    """Distance from center of storm"""
    track_data = track.get_track(cd.frameno)

    if track_data[0] is not None and track_data[1] is not None:
        return np.sqrt((cd.x - track_data[0])**2 + (cd.y - track_data[1])**2)
    else:
        return None


# ========================================================================
#  Water helper functions
# ========================================================================
def b(cd):
    return cd.aux[bathy_index, :, :]


def extract_eta(h, eta, DRY_TOL=1e-3):
    index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
    eta[index[0], index[1]] = np.nan
    return eta


def extract_velocity(h, hu, DRY_TOL=1e-8):
    u = np.zeros(hu.shape)
    index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
    u[index[0], index[1]] = hu[index[0], index[1]] / h[index[0], index[1]]
    return u


def eta(cd):
    return extract_eta(cd.q[0, :, :], cd.q[3, :, :])


def water_u(cd):
    return extract_velocity(cd.q[0, :, :], cd.q[1, :, :])


def water_v(cd):
    return extract_velocity(cd.q[0, :, :], cd.q[2, :, :])


def water_speed(current_data):
    u = water_u(current_data)
    v = water_v(current_data)

    return np.sqrt(u**2+v**2)


# ========================================================================
#  Plot items
# ========================================================================
def add_surface_elevation(plotaxes, plot_type='pcolor', bounds=None,
                          contours=None, shrink=1.0):
    """Add plotitem representing the sea surface."""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='surface', plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.surface_or_depth

        if bounds is not None:
            if bounds[0] == 0.0:
                plotitem.pcolor_cmap = plt.get_cmap('OrRd')
            else:
                plotitem.pcolor_cmap = surface_cmap
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Surface Height (m)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 0, 0, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='surface',
                                         plot_type='2d_contour')
        if bounds is None:
            plotitem.contour_levels = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]

        plotitem.plot_var = geoplot.surface_or_depth
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 0, 0, 0]
        plotitem.amr_contour_colors = 'k'

    elif plot_type == 'contourf':
        plotitem = plotaxes.new_plotitem(name='surface',
                                         plot_type='2d_contourf')
        plotitem.plot_var = geoplot.surface_or_depth
        if bounds is not None:
            contours = numpy.linspace(bounds[0], bounds[1], 11)
            plotitem.contour_levels = contours
            plotitem.fill_cmin = bounds[0]
            plotitem.fill_cmax = bounds[1]
        elif contours is not None:
            plotitem.contour_levels = contours
            plotitem.fill_cmin = min(contours)
            plotitem.fill_cmax = max(contours)

        plotitem.add_colorbar = True
        plotitem.fill_cmap = geoplot.tsunami_colormap
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Surface Height (m)"
        plotitem.fill_cmap = plt.get_cmap('OrRd')
        if any((value < 0 for value in plotitem.contour_levels)):
            plotitem.fill_cmap = surface_cmap

        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 0, 0, 0]
        plotitem.amr_contour_colors = 'k'


def add_speed(plotaxes, plot_type='pcolor', bounds=None,  contours=None,
              shrink=1.0):
    """Add plotitem representing speed of the water."""
    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_pcolor')
        plotitem.plot_var = water_speed
        # plotitem.plot_var = 1
        plotitem.pcolor_cmap = speed_cmap
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Current (m/s)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 0, 0, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_contour')
        if bounds is None:
            plotitem.contour_levels = [0.5, 1.5, 3, 4.5, 6.0]
        plotitem.kwargs = {'linewidths': 1}

        plotitem.plot_var = water_speed
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]
        plotitem.amr_contour_colors = 'k'

    elif plot_type == 'contourf':

        plotitem = plotaxes.new_plotitem(name='speed', plot_type='2d_contourf')

        plotitem.add_colorbar = True
        plotitem.colorbar_label = "Current (m/s)"
        plotitem.colorbar_shrink = shrink
        plotitem.fill_cmap = plt.get_cmap('PuBu')
        if bounds is not None:
            plotitem.contour_levels = numpy.linspace(bounds[0], bounds[1], 11)
            plotitem.fill_cmin = bounds[0]
            plotitem.fill_cmap = bounds[1]
        elif contours is not None:
            plotitem.contour_levels = contours
            plotitem.fill_cmin = min(contours)
            plotitem.fill_cmax = max(contours)

        # Modify the 'extends' plot attribute as we don't want this to extend
        # below 0
        plotitem.kwargs['extend'] = 'max'

        plotitem.plot_var = water_speed
        plotitem.amr_contour_show = [1] * 10
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]
        plotitem.amr_contour_colors = 'k'


def add_friction(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    """Add plotitem for the friction field"""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='friction',
                                         plot_type='2d_pcolor')
        plotitem.plot_var = friction
        plotitem.pcolor_cmap = friction_cmap
        plotitem.colorbar_shrink = shrink
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Manning's-$n$ Coefficient"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [0] * 10


def add_wind(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    """Add plotitem for the wind speed."""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name='wind', plot_type='2d_pcolor')
        plotitem.plot_var = wind_speed
        plotitem.pcolor_cmap = wind_cmap
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Wind Speed (m/s)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name='wind', plot_type='2d_contour')
        plotitem.plot_var = wind_speed
        plotitem.contour_nlevels = len(surge_data.wind_refine)
        plotitem.countour_min = surge_data.wind_refine[0]
        plotitem.patchedges_show = 1


def add_pressure(plotaxes, bounds=None, plot_type='pcolor', shrink=1.0):
    """Add plotitem for the pressure field."""

    if plot_type == 'pcolor' or plot_type == 'imshow':
        plotitem = plotaxes.new_plotitem(name="pressure",
                                         plot_type='2d_pcolor')
        plotitem.plot_var = pressure
        plotitem.colorbar_shrink = shrink
        plotitem.pcolor_cmap = pressure_cmap
        if bounds is not None:
            plotitem.pcolor_cmin = bounds[0]
            plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = True
        plotitem.colorbar_shrink = shrink
        plotitem.colorbar_label = "Pressure (mbar)"
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]
    elif plot_type == 'contour':
        pass


def add_land(plotaxes, plot_type='pcolor', bounds=[0, 50]):
    """Add plotitem for land"""

    if plot_type == 'pcolor':
        plotitem = plotaxes.new_plotitem(name='land', plot_type='2d_pcolor')
        plotitem.show = True
        plotitem.plot_var = geoplot.land
        plotitem.pcolor_cmap = land_cmap
        plotitem.pcolor_cmin = bounds[0]
        plotitem.pcolor_cmax = bounds[1]
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [0] * 10
        plotitem.amr_patchedges_show = [1, 1, 1, 1, 1, 0, 0]

    elif plot_type == 'contour':
        plotitem = plotaxes.new_plotitem(name="land", plot_type='2d_contour')
        plotitem.plot_var = geoplot.land
        plotitem.contour_nlevels = 40
        plotitem.contour_min = bounds[0]
        plotitem.contour_max = bounds[1]
        plotitem.amr_contour_colors = ['g']  # color on each level
        plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0


def add_bathy_contours(plotaxes, contour_levels=None, color='k'):
    """Add plotitem to plot contours of the topography"""

    plotitem = plotaxes.new_plotitem(name='bathy', plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    if contour_levels is None:
        contour_levels = [0.0]
    plotitem.contour_levels = contour_levels
    plotitem.amr_contour_colors = [color]
    plotitem.kwargs = {'linestyles': 'solid', 'linewidths': 2}
    plotitem.amr_contour_show = [1] * 10
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


def add_storm_radii(plotaxes, track, radii=[100e3], color='r'):
    """Add radii to plots based on storm position"""
    plotitem = plotaxes.new_plotitem(name="storm radius", 
                                     plot_type="2d_contour")
    plotitem.plot_var = lambda cd: storm_radius(cd, track)
    plotitem.contour_levels = radii
    plotitem.contour_colors = color


# ===== Storm related plotting =======
def sec2days(seconds):
    """Converst seconds to days."""
    return seconds / (60.0**2 * 24.0)


# def plot_track(t, x, y, wind_radius, wind_speed, Pc, name=None):
#     r"""Plot hurricane track given a storm.data file"""

#     if name is None:
#         name = ""
#     else:
#         name = " - %s" % name

#     colors = ['r', 'b']
#     divide = (np.max(Pc) + np.min(Pc)) / 2.0

#     fig = plt.figure(1)
#     axes = fig.add_subplot(111)
#     indices = Pc < divide
#     axes.scatter(x[indices], y[indices], color='r', marker='o')
#     indices = Pc >= divide
#     axes.scatter(x[indices], y[indices], color='b', marker='o')
#     axes.set_title("Track%s" % name)

#     fig = plt.figure(2, figsize=(24, 6))
#     axes = fig.add_subplot(131)
#     axes.plot(sec2days(t), wind_speed)
#     axes.set_title("Maximum Wind Speed%s" % name)

#     axes = fig.add_subplot(132)
#     axes.plot(sec2days(t), wind_radius)
#     axes.set_title("Maximum Wind Radius%s" % name)

#     axes = fig.add_subplot(133)
#     axes.plot(sec2days(t), Pc)
#     axes.plot(sec2days(t), np.ones(t.shape) * divide, 'k--')
#     axes.set_title("Central Pressure%s" % name)


"""Plot the track and optionally the intensity of the storm

Easily plot the track and intensity of a storm using a mapping package.

:Input:
    - *axes* (matplotlib.pyplot.axes) Axes to plot into.  Default is *None*
    - *intensity* (bool) Plot the intensity of storm along the track.
    Defaults to *False*.
    - *track_color* (str) String or specification of plotting color to use
    for the track if *intensity* is not being plotted.
    - *category_color* (dict) Dictionary containing mapping between
    category numerical value and colors.  Defaults to [-1, 5] -> ['gray',
    'black', 'violet', 'blue', 'yellow', 'orange', 'red']
    - *categorization* (str) Type of categorization, always use *"NHC"*
- *legend_loc* (str) Location of legend. Available options are 'best' (0), 'upper right' (1),
    'upper left'(2), 'lower left'(3), 'lower right'(4), 'right'(5), 'center left'(6), 'center right'(7), 
    'lower center'(8), 'upper center'(9), 'center'(10). Default is 'best'.

:Output:
    - (matplotlib.pyplot.axes) Axes object that was plotted into.
"""
# ========================================================================
#  Returns axes
#  Storm with category plotting function
# ========================================================================


# def add_track(Storm, axes, plot_package=None, category_color=None, legend_loc='best',
#               intensity=False, categorization="NHC", limits=None, track_color='red'):

#     if category_color is None:
#         category_color = {5: 'red',
#                           4: 'orange',
#                           3: 'yellow',
#                           2: 'blue',  # edit color
#                           1: 'violet',
#                           0: 'black',
#                           -1: 'gray'}
#     category = Storm.category(categorization=categorization)

#     # make it if intensity = true

#     # basic plotting
#     longitude = Storm.eye_location[:, 0]
#     latitude = Storm.eye_location[:, 1]
#     for i in range(len(longitude)):
#         if intensity:
#             color = category_color[category[i]]
#         else:
#             color = track_color
#         axes.plot(longitude[i:i + 2], latitude[i:i + 2], color=color)

#     axes.set_xlabel("Longitude")
#     axes.set_ylabel("Latitude")

#     categories_legend = []

#     if intensity and categorization == "NHC":
#         categories_legend = []
#         # plotitem = plotaxes.new_plotitem(name='category', plot_type='1d_plot')

#         if (-1 in category):
#             negativeone = mlines.Line2D(
#                 [], [], color=category_color[-1], marker='s', ls='', label="Tropical Depression")
#             categories_legend.append(negativeone)

#         if (0 in category):
#             zero = mlines.Line2D(
#                 [], [], color=category_color[0], marker='s', ls='', label="Tropical Storn")
#             categories_legend.append(zero)

#         if (1 in category):
#             one = mlines.Line2D([], [], color=category_color[1],
#                                 marker='s', ls='', label="Category 1")
#             categories_legend.append(one)

#         if (2 in category):
#             two = mlines.Line2D([], [], color=category_color[2],
#                                 marker='s', ls='', label="Category 2")
#             categories_legend.append(two)

#         if (3 in category):
#             three = mlines.Line2D(
#                 [], [], color=category_color[3], marker='s', ls='', label="Category 3")
#             categories_legend.append(three)

#         if (4 in category):
#             four = mlines.Line2D(
#                 [], [], color=category_color[4], marker='s', ls='', label="Category 4")
#             categories_legend.append(four)

#         if (5 in category):
#             five = mlines.Line2D(
#                 [], [], color=category_color[5], marker='s', ls='', label="Category 5")
#             categories_legend.append(five)

#         plt.legend(handles=categories_legend, loc=legend_loc)

#     # if bounds is not None:
#     #     plotitem.pcolor_cmin = bounds[0]
#     #     plotitem.pcolor_cmax = bounds[1]

#     return axes
