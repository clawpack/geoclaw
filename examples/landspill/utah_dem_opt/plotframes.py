#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the BSD 3-Clause license.

"""
Plot topography and flow depth
"""

import os
import numpy
import matplotlib
from matplotlib import pyplot
from clawpack.geoclaw import topotools
from clawpack import pyclaw


def get_max_depth(solution):
    """
    """

    max_depth = 0.

    for state in solution.states:
        max_temp = state.q[0, :, :].max()
        if max_temp > max_depth:
            max_depth = max_temp

    return max_depth


def get_max_AMR_level(solution):
    """
    Get the max AMR level in a solution object.
    """

    max_level = 1

    for state in solution.states:
        p = state.patch

        if p.level > max_level:
            max_level = p.level

    return max_level


def get_level_ncells_volumes(solution):
    """
    Get level-wise numbers of cells and fluid volumes.
    """

    max_level = get_max_AMR_level(solution)

    ncells = numpy.zeros(max_level, dtype=numpy.int)
    volumes = numpy.zeros(max_level, dtype=numpy.float64)

    for state in solution.states:
        p = state.patch
        level = p.level

        ncells[level-1] += p.num_cells_global[0] * p.num_cells_global[1]
        volumes[level-1] += (numpy.sum(state.q[0, :, :]) * p.delta[0] * p.delta[1])

    return ncells, volumes


def plot_at_axes(solution, ax, min_level=2, max_level=None,
                 shift=[0., 0.], vmin=0, vmax=None, dry_tol=1e-4,
                 cmap=pyplot.cm.viridis):
    """
    Plot fluid depth
    """

    if vmax is None:
        vmax = get_max_depth(solution)

    for state in solution.states:
        p = state.patch

        if p.level < min_level:
            continue

        if max_level is not None:
            if p.level > max_level:
                continue

        x, dx = numpy.linspace(p.lower_global[0], p.upper_global[0],
                               p.num_cells_global[0]+1, retstep=True)
        y, dy = numpy.linspace(p.lower_global[1], p.upper_global[1],
                               p.num_cells_global[1]+1, retstep=True)
        assert numpy.abs(dx-p.delta[0]) < 1e-6, "{} {}".format(dx, p.delta[0])
        assert numpy.abs(dy-p.delta[1]) < 1e-6, "{} {}".format(dy, p.delta[1])

        x -= shift[0]
        y -= shift[1]

        im = ax.pcolormesh(
            x, y, numpy.ma.masked_less(state.q[0, :, :], dry_tol).T,
            shading='flat', edgecolors='None',
            vmin=vmin, vmax=vmax, cmap=cmap)

    try:
        return im
    except UnboundLocalError:
        return None


topo_xbg = -12461650
topo_xed = topo_xbg + 4000
topo_ybg = 4984000
topo_yed = topo_ybg + 4000

source_x = topo_xbg + 2000
source_y = topo_ybg + 2000

domain_xbg = source_x - 400
domain_xed = source_x + 400
domain_ybg = source_y - 400
domain_yed = source_y + 400

raster_x_size = 500
raster_y_size = 350

x_crop_bg = source_x - raster_x_size / 2
x_crop_ed = x_crop_bg + raster_x_size
y_crop_bg = source_y - raster_y_size / 2
y_crop_ed = y_crop_bg + raster_y_size


# read topo file
topofile = topotools.Topography()
topofile.read("./utah_dem_topo_3.txt", topo_type=3)
topofile_crop = topofile.crop([x_crop_bg, x_crop_ed, y_crop_bg, y_crop_ed])

colormap = pyplot.cm.terrain
topo_min = 1284
topo_max = 1290
frame_bg = 0
frame_ed = 240

# read and plot solution with pyclaw
for frameno in range(frame_bg, frame_ed+1):

    if not os.path.isdir("./_plots"):
        os.mkdir("_plots")

    if os.path.isfile("./_plots/frame{:04d}.png".format(frameno)):
        print("Frame No. {} exists. Skip.".format(frameno))
        continue

    # a new figure
    fig = pyplot.figure(0, (13, 8), 90)

    # create an axes at 1, 3, 1
    ax_topo = fig.add_axes([0.1, 0.125, 0.65, 0.75])

    # light source
    ls = matplotlib.colors.LightSource(315, 45)

    # show topography in cropped region
    ax_topo.imshow(
        ls.shade(
            topofile_crop.Z,
            blend_mode='overlay',
            vert_exag=3,
            dx=1, dy=1, vmin=topo_min, vmax=topo_max,
            cmap=colormap),
        origin='lower')

    # set x and y ticks
    ax_topo.set_xticks(numpy.linspace(0, raster_x_size, 6))
    ax_topo.set_xticklabels(numpy.linspace(x_crop_bg, x_crop_ed, 6), rotation=-45, ha="left")
    ax_topo.set_yticks(numpy.linspace(0, raster_y_size, 8))
    ax_topo.set_yticklabels(numpy.linspace(y_crop_bg, y_crop_ed, 8))

    # x, y labels
    ax_topo.set_xlabel("x coordinates (m)")
    ax_topo.set_ylabel("y coordinates (m)")


    # plot colorbar in a new axes for topography
    cbarax = fig.add_axes([0.775, 0.125, 0.03, 0.75])
    im = ax_topo.imshow(topofile_crop.Z, cmap=colormap,
                        vmin=topo_min, vmax=topo_max, origin='lower')
    im.remove()
    cbar = pyplot.colorbar(im, cax=cbarax, ax=ax_topo)
    cbar.set_label("Elevation (m)")

    soln = pyclaw.Solution()

    aux = False
    if os.path.isfile("./_output/fort.a" + "{}".format(frameno).zfill(4)):
        aux = True

    soln.read(frameno, file_format="binary", read_aux=aux)

    ncells, volumes = get_level_ncells_volumes(soln)

    print("Plotting frame No. {}, T={} secs ({} mins)".format(
        frameno, soln.state.t, int(soln.state.t/60.)))

    print("\tFrame No. {}: N Cells: {}; Fluid volume: {}".format(
        frameno, ncells, volumes))

    im = plot_at_axes(soln, ax_topo, min_level=2,
                      shift=[x_crop_bg, y_crop_bg], dry_tol=1e-8)

    ax_topo.set_xlim(0, raster_x_size)
    ax_topo.set_ylim(0, raster_y_size)

    # plot colorbar in a new axes for depth
    cbarax = fig.add_axes([0.875, 0.125, 0.03, 0.75])
    if im is None:
        im = ax_topo.pcolormesh([0, raster_x_size], [0, raster_y_size], [[0]])
        im.remove()
        cbar = pyplot.colorbar(im, cax=cbarax, ax=ax_topo)
        cbar.ax.set_yticklabels([0]*len(cbar.ax.get_yticks()))
    else:
        cbar = pyplot.colorbar(im, cax=cbarax, ax=ax_topo)
    cbar.set_label("Depth (m)")

    # plot point source
    line = ax_topo.plot(source_x, source_y, 'r.', markersize=10)

    # figure title
    fig.suptitle("Topography and depth near breakage point, "
                 "T = {} (mins)".format(int(soln.state.t/60.)),
                 x=0.5, y=0.9, fontsize=16,
                 horizontalalignment="center",
                 verticalalignment="bottom")

    # label
    legend = ax_topo.legend(line, ["source"])

    # save image
    fig.savefig("./_plots/frame{:04d}.png".format(frameno), dpi=90)

    # clear
    pyplot.close(fig)
