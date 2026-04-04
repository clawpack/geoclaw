#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing GeoClaw data files."""


from pathlib import Path
import numpy as np
import pytest
import clawpack.geoclaw.data
import clawpack.geoclaw.fgmax_tools as fgmax_tools


def _read_text(path):
    """Read a text file for simple content checks in round-trip tests."""
    return Path(path).read_text()


@pytest.mark.python
def test_read_fgmax_data(tmp_path):
    r"""Test reading and writing of FGmaxData files."""
    data_file = tmp_path / "fgmax_grids.data"

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
    fg.tstart_max = 10.0
    fg.tend_max = 1.0e10
    fg.dt_check = 60.0
    fg.min_level_check = 3
    fg.arrival_tol = 1.0e-2
    fg.interp_method = 0
    fgmax_data.fgmax_grids.append(fg)

    fgmax_data.write(out_file=data_file)

    # Read data object
    read_fgmax_data = clawpack.geoclaw.data.FGmaxData()
    read_fgmax_data.read(data_file)

    assert read_fgmax_data.num_fgmax_val == fgmax_data.num_fgmax_val
    assert len(read_fgmax_data.fgmax_grids) == 1

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


# Additional FGmaxData round-trip test with multiple grids and point styles
@pytest.mark.python
@pytest.mark.xfail(reason="FGmaxData.read does not yet round-trip this multi-grid/point_style case correctly.")
def test_read_fgmax_data_multiple_grids(tmp_path):
    r"""Test FGmaxData round-trip with multiple grids and point styles."""
    data_file = tmp_path / "fgmax_grids_multi.data"

    fgmax_data = clawpack.geoclaw.data.FGmaxData()
    fgmax_data.num_fgmax_val = 1

    fg1 = fgmax_tools.FGmaxGrid()
    fg1.point_style = 2
    fg1.dx = 0.25
    fg1.x1 = -1.0
    fg1.x2 = 1.0
    fg1.y1 = -2.0
    fg1.y2 = 0.0
    fg1.tstart_max = 0.0
    fg1.tend_max = 100.0
    fg1.dt_check = 10.0
    fg1.min_level_check = 1
    fg1.arrival_tol = 1.0e-3
    fg1.interp_method = 0

    fg2 = fgmax_tools.FGmaxGrid()
    fg2.point_style = 1
    fg2.npts = 3
    fg2.xy_fname = "fgmax_points.txt"
    fg2.tstart_max = 5.0
    fg2.tend_max = 50.0
    fg2.dt_check = 5.0
    fg2.min_level_check = 2
    fg2.arrival_tol = 5.0e-3
    fg2.interp_method = 1

    fgmax_data.fgmax_grids.extend([fg1, fg2])
    fgmax_data.write(out_file=data_file)

    read_fgmax_data = clawpack.geoclaw.data.FGmaxData()
    read_fgmax_data.read(data_file)

    assert read_fgmax_data.num_fgmax_val == fgmax_data.num_fgmax_val
    assert len(read_fgmax_data.fgmax_grids) == 2

    rfg1, rfg2 = read_fgmax_data.fgmax_grids
    assert rfg1.point_style == fg1.point_style
    assert np.allclose(rfg1.dx, fg1.dx)
    assert np.allclose(rfg1.x1, fg1.x1)
    assert np.allclose(rfg1.x2, fg1.x2)
    assert np.allclose(rfg1.y1, fg1.y1)
    assert np.allclose(rfg1.y2, fg1.y2)

    assert rfg2.point_style == fg2.point_style
    assert rfg2.npts == fg2.npts
    assert rfg2.xy_fname == fg2.xy_fname
    assert np.allclose(rfg2.tstart_max, fg2.tstart_max)
    assert np.allclose(rfg2.tend_max, fg2.tend_max)
    assert np.allclose(rfg2.dt_check, fg2.dt_check)
    assert np.allclose(rfg2.arrival_tol, fg2.arrival_tol)
    assert np.allclose(rfg2.interp_method, fg2.interp_method)


# DTopoData round-trip test
@pytest.mark.python
@pytest.mark.xfail(reason="DTopoData.read currently fails to parse this written dtopo.data round-trip case.")
def test_dtopo_data_roundtrip(tmp_path):
    r"""Test reading and writing of DTopoData files."""
    data_file = tmp_path / "dtopo.data"

    dtopo_data = clawpack.geoclaw.data.DTopoData()
    dtopo_data.dt_max_dtopo = 2.5
    dtopo_data.dtopofiles = [
        [1, 2, 3, "dtopo_one.tt3"],
        [3, 4, 1, "dtopo_two.tt1"],
    ]

    dtopo_data.write(out_file=data_file)

    read_dtopo_data = clawpack.geoclaw.data.DTopoData()
    read_dtopo_data.read(data_file)

    assert np.allclose(read_dtopo_data.dt_max_dtopo, dtopo_data.dt_max_dtopo)
    assert read_dtopo_data.dtopofiles == dtopo_data.dtopofiles

    text = _read_text(data_file)
    assert "dtopo_one.tt3" in text
    assert "dtopo_two.tt1" in text


# SurgeData round-trip test
@pytest.mark.python
def test_surge_data_roundtrip(tmp_path):
    r"""Test reading and writing of SurgeData files."""
    data_file = tmp_path / "surge.data"

    surge_data = clawpack.geoclaw.data.SurgeData()
    surge_data.wind_forcing = True
    surge_data.drag_law = 2
    surge_data.pressure_forcing = True
    surge_data.wind_index = 6
    surge_data.pressure_index = 7
    surge_data.display_landfall_time = True
    surge_data.wind_refine = [20.0, 40.0, 60.0]
    surge_data.R_refine = [60.0e3, 40.0e3, 20.0e3]
    surge_data.storm_specification_type = "data"
    surge_data.storm_file = "synthetic.storm"

    surge_data.write(out_file=data_file)

    read_surge_data = clawpack.geoclaw.data.SurgeData()
    read_surge_data.read(data_file)

    assert read_surge_data.wind_forcing == surge_data.wind_forcing
    assert read_surge_data.drag_law == surge_data.drag_law
    assert read_surge_data.pressure_forcing == surge_data.pressure_forcing
    assert read_surge_data.wind_index == surge_data.wind_index
    assert read_surge_data.pressure_index == surge_data.pressure_index
    assert read_surge_data.display_landfall_time == surge_data.display_landfall_time
    assert np.allclose(read_surge_data.wind_refine, surge_data.wind_refine)
    assert np.allclose(read_surge_data.wind_refine, surge_data.wind_refine)
    assert np.allclose(read_surge_data.R_refine, surge_data.R_refine)
    expected_spec = clawpack.geoclaw.data.SurgeData.storm_spec_dict_mapping["data"]
    assert read_surge_data.storm_specification_type == expected_spec
    assert read_surge_data.storm_file == surge_data.storm_file

    text = _read_text(data_file)
    assert "synthetic.storm" in text


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
