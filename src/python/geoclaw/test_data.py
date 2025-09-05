#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing GeoClaw data files"""

from pathlib import Path
import tempfile
import shutil

import numpy as np

import clawpack.clawutil.test as test
import clawpack.geoclaw.data
import clawpack.geoclaw.fgmax_tools as fgmax_tools


def test_read_FGmaxData(tmp_path):
    r"""Test readinga and writing of FGmaxData files"""

    try:
        data_path = Path(tmp_path) / "fgmax_grids.data"

        # Test data object
        fgmax_data = clawpack.geoclaw.data.FGmaxData()
        fgmax_data.num_fgmax_val = 2

        # Test grid data
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2
        fg.dx = 2.0 / (3.0 * 4.0)
        fg.x1 = -120.0 + fg.dx / 2.0
        fg.x2 = -60.0 - fg.dx / 2.0
        fg.y1 = -60.0 + fg.dx / 2.0
        fg.y2 = 0.0 - fg.dx / 2.0
        fg.tstart_max =  10.0
        fg.tend_max = 1.e10
        fg.dt_check = 60.0
        fg.min_level_check = 3
        fg.arrival_tol = 1.e-2
        fg.interp_method = 0
        fgmax_data.fgmax_grids.append(fg)
    
        fgmax_data.write(out_file=data_path)

        # Read data object
        read_fgmax_data = clawpack.geoclaw.data.FGmaxData()
        read_fgmax_data.read(data_path)

        # Tests
        tfg = read_fgmax_data.fgmax_grids[0]
        assert np.allclose(fg.x1, tfg.x1)
        assert np.allclose(fg.x2, tfg.x2)
        assert np.allclose(fg.y1, tfg.y1)
        assert np.allclose(fg.y2, tfg.y2)
        assert np.allclose(fg.tstart_max, tfg.tstart_max)
        assert np.allclose(fg.tend_max, tfg.tend_max)
        assert np.allclose(fg.dt_check, tfg.dt_check)
        assert np.allclose(fg.min_level_check, tfg.min_level_check)
        assert np.allclose(fg.arrival_tol, tfg.arrival_tol)
        assert np.allclose(fg.interp_method, tfg.interp_method)


    except Exception as e:
        print("Error")
        path = Path() / "test_read_FGmaxData_output"
        shutil.rmtree(path, ignore_errors=True)
        shutil.copytree(tmp_path, path)
        raise e
