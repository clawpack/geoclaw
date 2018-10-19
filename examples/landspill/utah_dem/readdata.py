#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
#
# Distributed under terms of the MIT license.

"""
Read output ASCII data
"""
import numpy
import scipy.interpolate
import matplotlib.pyplot as pyplot

class SolutionInfo:
    """
    Class to hold information of a solution file
    """

    def __init__(self, filename):
        """
        initialize this instance wiht a file, fort.tXXXX
        """

        raw_data = numpy.genfromtxt(filename, dtype=None, encoding='utf8')

        for row in raw_data:
            self.__dict__[row[1]] = row[0]

        del raw_data

class SolutionPatch:
    """
    Solution pathc
    """

    def __init__(self):
        """
        initialization
        """

        self.grid_number = None
        self.AMR_level = None
        self.mx = None
        self.my = None
        self.xlow = None
        self.ylow = None
        self.dx = None
        self.dy = None
        self.data = None
        self.interp = None

        self.xhigh = None
        self.yhigh = None

        self.valid = False

    def __str__(self):
        """
        __str__
        """
        return "grid_number: {}\n".format(self.grid_number) + \
            "AMR_level: {}\n".format(self.AMR_level) + \
            "mx: {}\n".format(self.mx) + \
            "my: {}\n".format(self.my) + \
            "xlow: {}\n".format(self.xlow) + \
            "ylow: {}\n".format(self.ylow) + \
            "xhigh: {}\n".format(self.xhigh) + \
            "yhigh: {}\n".format(self.yhigh) + \
            "dx: {}\n".format(self.dx) + \
            "dy: {}\n".format(self.dy)

    def complete_info(self):
        """
        Complete the information of this patch
        """
        self.xhigh = self.xlow + self.mx * self.dx
        self.yhigh = self.ylow + self.my * self.dy

    def check(self):
        """
        check status
        """
        assert self.grid_number is not None
        assert self.AMR_level is not None
        assert self.mx is not None
        assert self.my is not None
        assert self.xlow is not None
        assert self.ylow is not None
        assert numpy.abs(self.xhigh-self.xlow-self.dx*self.mx) < 1e-6
        assert numpy.abs(self.yhigh-self.ylow-self.dy*self.my) < 1e-6
        assert self.dx is not None
        assert self.dy is not None

        if self.valid:
            assert self.data is not None
            assert self.data.shape == (self.my, self.mx, 4)
        else:
            assert self.data is None

    def _create_interp(self):
        """
        Create scipy.interpolate.RectBivariateSpline instance
        """

        self.check()

        if self.interp is not None:
            return

        x, dx = numpy.linspace(self.xlow+self.dx/2., self.xhigh-self.dx/2.,
                               self.mx, retstep=True)
        y, dy = numpy.linspace(self.ylow+self.dy/2., self.yhigh-self.dy/2.,
                               self.my, retstep=True)
        assert numpy.abs(dx-self.dx) < 1e-6, "{} {}".format(dx, self.dx)
        assert numpy.abs(dy-self.dy) < 1e-6, "{} {}".format(dy, self.dy)

        self.interp = numpy.empty(4, dtype=scipy.interpolate.RectBivariateSpline)

        for i in range(4):
            self.interp[i] = scipy.interpolate.RectBivariateSpline(
                x, y, self.data[:, :, i].T,
                [self.xlow, self.xhigh, self.ylow, self.yhigh])

    def interpolate(self, x, y, shift=[0., 0.], pos=0):
        """
        Interpolate with 3rd spline
        """
        self._create_interp()

        return self.interp[pos](x-shift[0], y-shift[1]).T

class Solution(SolutionInfo):
    """
    Class to hold solution
    """

    def __init__(self, outdir, frameno, min_lvl=1, max_lvl=None):
        """
        initialization
        """

        super().__init__(outdir+"/fort.t{:04d}".format(frameno))

        self.interp = numpy.empty(4, dtype=scipy.interpolate.interp2d)
        self.interp[:] = None
        self.max_dpeth = -99999.
        self.min_lvl = min_lvl

        with open(outdir+"/fort.q{:04d}".format(frameno), encoding='utf8') as f:
            lines = f.readlines()

        int_keys = ["grid_number", "AMR_level", "mx", "my"]
        float_keys = ["xlow", "ylow", "dx", "dy"]

        patch_count = 0
        local_header_count = 0
        local_data_count = 0
        global_line_count = 0
        self.patches = numpy.empty(0, dtype=SolutionPatch)

        while True:

            # split the current line
            line = lines[global_line_count].strip().split()

            # this line is part of the header of a patch
            if len(line) == 2:

                # make sure it's legal header info
                assert (line[1] in int_keys) or (line[1] in float_keys)

                # this means it's the beginning line of a new patch
                if local_header_count == 0:
                    self.patches = numpy.append(self.patches, SolutionPatch())

                if line[1] in int_keys:
                    self.patches[-1].__dict__[line[1]] = numpy.int(line[0])

                if line[1] in float_keys:
                    self.patches[-1].__dict__[line[1]] = numpy.float(line[0])

                local_header_count += 1

                # if all header info has been read for this patch
                if local_header_count == 8:

                    # check if this pactch is in the desired range
                    if (self.patches[-1].AMR_level >= min_lvl):
                        if (max_lvl is None) or (
                                self.patches[-1].AMR_level <= max_lvl):
                            self.patches[-1].valid = True

                    if self.patches[-1].valid:
                        self.patches[-1].data = numpy.zeros(
                            (self.patches[-1].my, self.patches[-1].mx, 4),
                            dtype=numpy.float64)

                    self.patches[-1].complete_info()
                    self.patches[-1].check()
                    local_header_count = 0

                    # examinate the size of the data in the previous patch
                    try:
                        assert local_data_count == \
                            self.patches[-2].mx * self.patches[-2].my
                    except IndexError:
                        pass

                    # reset local_data_count
                    local_data_count = 0

            # data of the current pathc
            elif len(line) == 4:

                if self.patches[-1].valid:
                    i = local_data_count % self.patches[-1].mx
                    j = int(local_data_count / self.patches[-1].mx)
                    self.patches[-1].data[j, i] = [numpy.float64(v) for v in line]
                    if self.patches[-1].data[j, i, 0] > self.max_dpeth:
                        self.max_dpeth = self.patches[-1].data[j, i, 0]

                local_data_count += 1

            # empty line
            elif len(line) == 0:
                pass
            else:
                raise Error

            global_line_count += 1

            # break if all lines has looped over
            if global_line_count == len(lines):
                break

        # calculate max lvl
        self.max_lvl = 0
        for p in self.patches:
            if p.valid and p.AMR_level > self.max_lvl:
                self.max_lvl = p.AMR_level


    def plot_at_axes(self, ax, min_level=2, max_level=None,
                     shift=[0., 0.], vmin=0, vmax=None, dry_tol=1e-4,
                     cmap=pyplot.cm.viridis):
        """
        Plot fluid depth
        """
        for p in self.patches:

            if p.AMR_level < min_level:
                continue

            if max_level is not None:
                if p.AMR_level > max_level:
                    continue

            x, dx = numpy.linspace(p.xlow, p.xlow+p.mx*p.dx, p.mx+1, retstep=True)
            y, dy = numpy.linspace(p.ylow, p.ylow+p.my*p.dy, p.my+1, retstep=True)
            assert numpy.abs(dx-p.dx) < 1e-6, "{} {}".format(dx, p.dx)
            assert numpy.abs(dy-p.dy) < 1e-6, "{} {}".format(dy, p.dy)

            x -= shift[0]
            y -= shift[1]

            if vmax is None:
                vmax = self.max_dpeth

            im = ax.pcolormesh(
                x, y, numpy.ma.masked_less(p.data[:, :, 0], dry_tol),
                shading='flat', edgecolors='None',
                vmin=vmin, vmax=vmax, cmap=cmap)

        try:
            return im
        except UnboundLocalError:
            return None

    def interpolate(self, x, y, shift=[0., 0.], pos=0, level=2,
                    filter_less=1e-7, nodatavalue=-9999.):
        """
        interpolation
        """

        values = numpy.zeros((y.size, x.size), dtype=numpy.float64)

        for p in self.patches:

            if p.AMR_level != level:
                continue

            xid = numpy.where((x>=p.xlow)&(x<=p.xhigh))[0]
            yid = numpy.where((y>=p.ylow)&(y<=p.yhigh))[0]

            if xid.size and yid.size:
                values[yid[:, None], xid[None, :]] = p.interpolate(x[xid], y[yid])

        values[values<filter_less] = nodatavalue

        return values
