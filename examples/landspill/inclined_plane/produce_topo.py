#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the BSD 3-Clause license.

"""
Produce an inclined plane with degree 2.5
"""
import numpy

if __name__ == "__main__":
    xlow = -2.0
    ylow = -2.0
    Lx = 4.0
    Ly = 4.0
    cellsize = 0.005
    mx = int(Lx / cellsize + 0.5)
    my = int(Ly / cellsize + 0.5)

    x = numpy.linspace(xlow+cellsize/2., Lx+xlow-cellsize/2., mx)
    topo = numpy.zeros((my, mx), dtype=numpy.float64)

    for i in range(mx):
        topo[:, i] = (Lx + xlow - x[i]) * numpy.sin(2.5/180*numpy.pi)

    headers = "{}\t\t\tmx\n".format(mx) + \
        "{}\t\t\tmy\n".format(my) + \
        "{}\t\txlower\n".format(xlow) + \
        "{}\t\tylower\n".format(ylow) + \
        "{}\t\tcellsize\n".format(cellsize) + \
        "-9999\t\tnodatavalue\n"

    with open("inclined_plane_2.5.txt", "w") as f:
        f.write(headers)

        for j in range(my):
            topo[j, :].tofile(f, " ")
            f.write("\n")
