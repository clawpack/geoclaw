#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the BSD 3-Clause license.

"""
Write to NetCDF4 file with CF convention
"""

import os
import numpy
import scipy.interpolate
import netCDF4
import time as pytime
from osgeo import osr
from clawpack import pyclaw


def get_state_interpolator(state):
    """
    Get the scipy interpolation function of each patch.
    """
    p = state.patch

    x, dx = numpy.linspace(p.lower_global[0]+p.delta[0]/2.,
                           p.upper_global[0]-p.delta[0]/2.,
                           p.num_cells_global[0], retstep=True)
    y, dy = numpy.linspace(p.lower_global[1]+p.delta[1]/2.,
                           p.upper_global[1]-p.delta[1]/2.,
                           p.num_cells_global[1], retstep=True)
    assert numpy.abs(dx-p.delta[0]) < 1e-6, "{} {}".format(dx, p.delta[0])
    assert numpy.abs(dy-p.delta[1]) < 1e-6, "{} {}".format(dy, p.delta[1])

    interp = numpy.empty(4, dtype=scipy.interpolate.RectBivariateSpline)

    for i in range(4):
        interp[i] = scipy.interpolate.RectBivariateSpline(
            x, y, state.q[i, :, :],
            [p.lower_global[0], p.upper_global[0],
             p.lower_global[1], p.upper_global[1]])

    return interp


def interpolate(solution, x, y, shift=[0., 0.], pos=0, level=2,
                filter_less=1e-7, nodatavalue=-9999.):
    """
    Do the interpolation.
    """

    values = numpy.zeros((y.size, x.size), dtype=numpy.float64)

    for state in solution.states:
        p = state.patch

        if p.level != level:
            continue

        xid = numpy.where((x>=p.lower_global[0])&(x<=p.upper_global[0]))[0]
        yid = numpy.where((y>=p.lower_global[1])&(y<=p.upper_global[1]))[0]

        interpolator = get_state_interpolator(state)

        if xid.size and yid.size:
            values[yid[:, None], xid[None, :]] = \
                interpolator[0](x[xid]-shift[0], y[yid]-shift[1]).T

    values[values<filter_less] = nodatavalue

    return values


# some simulation parameters
source_x = -12459650
source_y = 4986000
out_Lx = 300
out_Ly = 175
out_dx = 1.
out_dy = 1.
data_dir = os.getcwd()

x_raw = numpy.arange(source_x-out_Lx+out_dx/2., source_x+out_Lx, out_dx)
y_raw = numpy.arange(source_y-out_Ly+out_dy/2., source_y+out_Ly, out_dy)

# create a NC file and root group
rootgrp = netCDF4.Dataset(
    "{}.nc".format(os.path.basename(os.getcwd())), mode="w", format="NETCDF4")

# set up dimension
ntime = rootgrp.createDimension("time", None)
nx = rootgrp.createDimension("x", x_raw.size)
ny = rootgrp.createDimension("y", y_raw.size)

# create variables
times = rootgrp.createVariable("time", numpy.float64, ("time",))
x = rootgrp.createVariable("x", numpy.float64, ("x",))
y = rootgrp.createVariable("y", numpy.float64, ("y",))
depth = rootgrp.createVariable(
    "depth", numpy.float64, ("time", "y", "x"),
    fill_value=-9999., zlib=True, complevel=9)
mercator = rootgrp.createVariable("mercator", 'S1')

# global attributes
rootgrp.title = os.path.basename(os.getcwd())
rootgrp.institution = "N/A"
rootgrp.source = "N/A"
rootgrp.history = "Created " + pytime.ctime(pytime.time())
rootgrp.reference = ""
rootgrp.comment = ""
rootgrp.Conventions = "CF-1.7"

# variable attributes: time
times.units = "sec"
times.axis = "T"
times.long_name = "Simulation times"

# variable attributes: x
x.units = "m"
x.long_name = "X-coordinate in EPSG:3857 WGS 84"
x.standard_name = "projection_x_coordinate"

# variable attributes: y
y.units = "m"
y.long_name = "Y-coordinate in EPSG:3857 WGS 84"
y.standard_name = "projection_y_coordinate"

# variable attributes: topo
depth.units = "m"
depth.long_name = "Oil Depth"
depth.grid_mapping = "mercator"

# variable attributes: mercator
srs = osr.SpatialReference()
srs.ImportFromEPSG(3857)
mercator.grid_mapping_name = "mercator"
mercator.long_name = "CRS definition"
mercator.spatial_ref = srs.ExportToWkt()
mercator.GeoTransform = "{} {} 0 {} 0 {}".format(
    source_x-out_Lx, out_dx, source_y+out_Ly, -out_dy)

# variable values
x[:] = x_raw
y[:] = y_raw

for frameno in range(0, 241):
    print("Frame No.", frameno)
    soln = pyclaw.Solution()

    aux = False
    if os.path.isfile("./_output/fort.a" + "{}".format(frameno).zfill(4)):
        aux = True

    soln.read(frameno, file_format="binary", read_aux=aux)
    times[frameno] = soln.state.t
    depth[frameno, :, :] = interpolate(soln, x_raw, y_raw)

rootgrp.close()
