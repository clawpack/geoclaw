#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the BSD 3-Clause license.

"""
Read output ASCII data
"""
import numpy
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from clawpack import pyclaw


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


if __name__ == "__main__":
    import os

    T = [32, 59, 122, 271, 486, 727]

    for i, t in enumerate(T):

        print("T =", t)

        aux = False
        if os.path.isfile("./_output/fort.a" + "{}".format(i+1).zfill(4)):
            aux = True
        soln = pyclaw.Solution(i+1, file_format="binary", read_aux=aux)

        max_level = get_max_AMR_level(soln)
        print("\tValid Max AMR level: ", max_level)

        ncells, volumes = get_level_ncells_volumes(soln)
        print("\tNumber of cells: ", ncells)
        print("\tTotal volume at each AMR level: ", volumes)

        fig = pyplot.figure()
        ax = fig.gca()
        im = plot_at_axes(soln, ax, max_level, vmax=0.005, cmap=pyplot.cm.jet)
        ax.set_axisbelow(True)
        ax.set_aspect('equal')
        ax.set_title("Silicon Oil Depth at T = {} sec".format(t))
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")

        lister = numpy.loadtxt("./lister_1992/T={}.csv".format(t),
                               delimiter=',', skiprows=1)
        pyplot.plot(lister[:, 0], lister[:, 1], 'ro',
                    markersize=5, label="Experimental data (Lister, 1992)")

        pyplot.xlim(-0.2, 1.0)
        pyplot.ylim(-0.3, 0.3)
        pyplot.grid(True, 'both', 'both')
        pyplot.legend(loc=0)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.02)
        cbar = pyplot.colorbar(im, cax=cax)
        cbar.set_label("Depth (m)")
        cbar.ax.tick_params(labelsize=8)

        if not os.path.isdir("./_plots"):
            os.mkdir("_plots")

        pyplot.savefig("./_plots/T=" + "{}".format(t).zfill(3) + "s.png")
        #pyplot.show()
