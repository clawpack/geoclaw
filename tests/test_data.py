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
def test_dtopo_data_roundtrip(tmp_path):
    r"""Test reading and writing of DTopoData files."""
    import clawpack.geoclaw.dtopotools as dtopotools

    data_file = tmp_path / "dtopo.data"

    d1 = dtopotools.DTopography()
    d1.path = "dtopo_one.tt3"
    d1.dtopo_type = 3
    d1.x_shift = 1.5
    d1.y_shift = -2.0
    d1.z_shift = -0.25
    d1.negate_z = True

    dtopo_data = clawpack.geoclaw.data.DTopoData()
    dtopo_data.dt_max_dtopo = 2.5
    dtopo_data.dtopofiles = [
        d1,
        [1, "dtopo_two.tt1"],   # legacy [dtopo_type, path] entry
    ]

    dtopo_data.write(out_file=data_file)

    read_dtopo_data = clawpack.geoclaw.data.DTopoData()
    read_dtopo_data.read(data_file)

    assert np.allclose(read_dtopo_data.dt_max_dtopo, dtopo_data.dt_max_dtopo)
    assert len(read_dtopo_data.dtopofiles) == 2

    r1, r2 = read_dtopo_data.dtopofiles
    assert r1.path.endswith("dtopo_one.tt3")
    assert r1.dtopo_type == 3
    assert np.allclose(r1.x_shift, 1.5)
    assert np.allclose(r1.y_shift, -2.0)
    assert np.allclose(r1.z_shift, -0.25)
    assert r1.negate_z is True
    assert r1.crop_extent is None
    assert r2.path.endswith("dtopo_two.tt1")
    assert r2.dtopo_type == 1
    assert r2.negate_z is False

    text = _read_text(data_file)
    assert "dtopo_one.tt3" in text
    assert "dtopo_two.tt1" in text


@pytest.mark.python
def test_dtopo_data_datum_mismatch_warns(tmp_path):
    r"""DTopoData.write warns when dtopo files carry different datums."""
    import clawpack.geoclaw.dtopotools as dtopotools

    a = dtopotools.DTopography()
    a.path = "a.tt3"
    a.dtopo_type = 3
    a.datum = "NAVD88"
    b = dtopotools.DTopography()
    b.path = "b.tt3"
    b.dtopo_type = 3
    b.datum = "MSL"

    dtopo_data = clawpack.geoclaw.data.DTopoData()
    dtopo_data.dtopofiles = [a, b]
    with pytest.warns(UserWarning, match="mismatched vertical datums"):
        dtopo_data.write(out_file=tmp_path / "dtopo.data")


@pytest.mark.python
def test_dtopo_data_no_datum_warning_when_consistent(tmp_path, recwarn):
    r"""No datum warning when dtopo datums agree or only one is set."""
    import clawpack.geoclaw.dtopotools as dtopotools

    a = dtopotools.DTopography()
    a.path = "a.tt3"
    a.dtopo_type = 3
    a.datum = "MSL"
    b = dtopotools.DTopography()
    b.path = "b.tt3"
    b.dtopo_type = 3            # datum stays None -> single distinct datum

    dtopo_data = clawpack.geoclaw.data.DTopoData()
    dtopo_data.dtopofiles = [a, b]
    dtopo_data.write(out_file=tmp_path / "dtopo.data")
    assert not any("mismatched vertical datums" in str(w.message)
                   for w in recwarn.list)


@pytest.mark.python
@pytest.mark.netcdf
def test_dtopo_data_netcdf_descriptor(tmp_path):
    r"""DTopoData.write() emits the descriptor block for type-4 entries."""
    pytest.importorskip("netCDF4")
    import clawpack.geoclaw.dtopotools as dtopotools

    # Build and write a small NetCDF dtopo file
    x = np.linspace(0.0, 2.0, 5)
    y = np.linspace(0.0, 1.0, 6)
    X, Y = np.meshgrid(x, y)
    src = dtopotools.DTopography()
    src.x, src.y, src.X, src.Y = x, y, X, Y
    src.times = [0.0, 0.5, 1.0]
    src.dZ = np.zeros((3,) + X.shape)
    nc_path = tmp_path / "dt.nc"
    src.write(nc_path, dtopo_type=4)

    entry = dtopotools.DTopography()
    entry.path = str(nc_path)
    entry.dtopo_type = 4

    dtopo_data = clawpack.geoclaw.data.DTopoData()
    dtopo_data.dtopofiles = [entry]
    data_file = tmp_path / "dtopo.data"
    dtopo_data.write(out_file=data_file)

    text = _read_text(data_file)
    assert "var_name       = dz" in text
    assert "x_name         = lon" in text
    assert "y_name         = lat" in text
    assert "t0             = 0.0" in text
    assert "dt             = 0.5" in text

    # The reader still parses the file (descriptor lines are skipped).
    read_back = clawpack.geoclaw.data.DTopoData()
    read_back.read(data_file)
    assert len(read_back.dtopofiles) == 1
    assert read_back.dtopofiles[0].dtopo_type == 4


@pytest.mark.python
def test_dtopo_data_unsupported_preprocessing(tmp_path):
    r"""Unsupported dtopo preprocessing attributes fail loudly at write."""
    import clawpack.geoclaw.dtopotools as dtopotools

    d = dtopotools.DTopography()
    d.path = "dtopo.tt3"
    d.dtopo_type = 3
    d.crop_extent = [0.0, 1.0, 0.0, 1.0]

    dtopo_data = clawpack.geoclaw.data.DTopoData()
    dtopo_data.dtopofiles = [d]

    with pytest.raises(NotImplementedError, match="crop_extent"):
        dtopo_data.write(out_file=tmp_path / "dtopo.data")


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
    surge_data.storm_time_scale = 2.5
    surge_data.t_ramp_on = 3600.0
    surge_data.t_ramp_off = 1800.0
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
    assert np.allclose(read_surge_data.storm_time_scale, surge_data.storm_time_scale)
    assert np.allclose(read_surge_data.t_ramp_on, surge_data.t_ramp_on)
    assert np.allclose(read_surge_data.t_ramp_off, surge_data.t_ramp_off)
    assert np.allclose(read_surge_data.wind_refine, surge_data.wind_refine)
    assert np.allclose(read_surge_data.wind_refine, surge_data.wind_refine)
    assert np.allclose(read_surge_data.R_refine, surge_data.R_refine)
    expected_spec = clawpack.geoclaw.data.SurgeData.storm_spec_dict_mapping["data"]
    assert read_surge_data.storm_specification_type == expected_spec
    assert read_surge_data.storm_file == surge_data.storm_file
    # Legacy "data" resolves to the gridded family on the new wire.
    assert read_surge_data.storm_family == "gridded"
    assert read_surge_data.storm_subtype == "gridded"

    text = _read_text(data_file)
    assert "synthetic.storm" in text


@pytest.mark.python
def test_surge_forcing_family_subtype(tmp_path):
    r"""family/subtype API, legacy aliases, and the fixed lookup bugs."""
    data = clawpack.geoclaw.data

    # Registry / alias resolution, including the previously-broken 'DeMaria'
    # case-mismatch and the hyphen/underscore + holland08/10 spellings.
    assert data.resolve_forcing_subtype("holland80") == ("parametric", "holland80", 1)
    assert data.resolve_forcing_subtype("DeMaria") == ("parametric", "demaria", 7)
    assert data.resolve_forcing_subtype("modified-rankine") == \
        ("parametric", "modified_rankine", 6)
    assert data.resolve_forcing_subtype("holland08") == ("parametric", "holland2008", 8)
    assert data.resolve_forcing_subtype("data") == ("gridded", "gridded", -1)
    assert data.resolve_forcing_subtype(None) == ("none", "none", 0)
    assert data.resolve_forcing_subtype(-1) == ("gridded", "gridded", -1)

    # New family/subtype API round-trips and populates the legacy selector.
    surge_data = data.SurgeData()
    surge_data.storm_family = "parametric"
    surge_data.storm_subtype = "holland2010"
    surge_data.storm_file = "b.storm"
    out = tmp_path / "surge.data"
    surge_data.write(out_file=out)

    read_back = data.SurgeData()
    read_back.read(out)
    assert read_back.storm_family == "parametric"
    assert read_back.storm_subtype == "holland2010"
    assert read_back.storm_specification_type == 2
    assert read_back.storm_file == "b.storm"

    # Explicit tokens are on the wire; the raw integer selector is not.
    text = _read_text(out)
    assert "storm_family" in text and "storm_subtype" in text
    assert "'parametric'" in text and "'holland2010'" in text

    # An inconsistent family/subtype pair is rejected at write time.
    bad = data.SurgeData()
    bad.storm_family = "gridded"
    bad.storm_subtype = "holland80"
    bad.storm_file = "x"
    with pytest.raises(ValueError, match="inconsistent"):
        bad.write(out_file=tmp_path / "bad.data")


# ---------------------------------------------------------------------------
# 1D topo.data / dtopo.data layout
#
# src/1d_classic shares these data classes with the 2D code but not the Fortran
# readers, and it was never updated for the per-file preprocessing block.  Its
# read_topo_settings expects
#
#     topo_missing / test_topography / ntopofiles / override_order / '<path>'
#
# and read_dtopo_settings reads the dtopo path and type with a single
# list-directed statement, so the path line must not carry a trailing comment.
# These pin both layouts; getting either wrong makes every 1d_classic example
# abort at startup.
# ---------------------------------------------------------------------------


def _payload(path):
    """Data-file lines with the generated comment header and blanks removed."""
    return [line for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]


def _field_names(path):
    """Ordered field labels of a .data payload, one per line.

    ``data_write`` emits ``<value>  =: <name>  # <description>``, while the
    path/type lines are written by hand as ``<value>   # <name>``.  Taking the
    ``=:`` label when present and the first comment token otherwise gives the
    file's layout as a list of names -- so a layout assertion can say *which*
    line is missing or misplaced rather than only that the count is wrong.
    """
    names = []
    for line in _payload(path):
        if "=:" in line:
            names.append(line.split("=:", 1)[1].split()[0])
        elif "#" in line:
            names.append(line.split("#", 1)[1].split()[0])
        else:
            names.append(line.strip())
    return names


@pytest.mark.python
def test_topo_data_1d_uses_legacy_layout(tmp_path):
    r"""1D topo.data carries override_order and no preprocessing block."""
    import clawpack.geoclaw.topotools as topotools

    topo = topotools.Topography()
    topo.path = "celledges.data"
    topo.topo_type = 1

    topo_data = clawpack.geoclaw.data.TopographyData(num_dim=1)
    topo_data.topofiles = [topo]
    out = tmp_path / "topo.data"
    topo_data.write(out_file=out)

    # Assert the layout by name, in order.  read_topo_settings reads these
    # positionally, so a dropped or reordered line is a startup abort -- and
    # naming them makes the failure say which one went missing.
    assert _field_names(out) == ["topo_missing", "test_topography",
                                 "ntopofiles", "override_order",
                                 "topo_path", "topo_type"]

    lines = _payload(out)
    # Unlike dtopo below, the topo path line may keep its trailing comment:
    # read_topo_settings reads the path as a *single* list-directed item, so
    # the read is satisfied before reaching the comment.
    assert lines[4].startswith("'") and "celledges.data'" in lines[4]
    text = Path(out).read_text()
    for attr in ("crop_extent", "coarsen", "buffer", "align", "negate_z"):
        assert attr not in text


@pytest.mark.python
def test_topo_data_1d_override_order_written_once_for_multiple_files(tmp_path):
    r"""override_order precedes the file entries and appears exactly once.

    read_topo_settings reads the logical *once*, before any path, so it must
    not migrate into the per-file loop.  With a single topo file a per-file
    write would be indistinguishable from the correct one; two files pin it.
    """
    import clawpack.geoclaw.topotools as topotools

    topos = []
    for name in ("coarse.tt3", "fine.tt3"):
        topo = topotools.Topography()
        topo.path = name
        topo.topo_type = 3
        topos.append(topo)

    topo_data = clawpack.geoclaw.data.TopographyData(num_dim=1)
    topo_data.topofiles = topos
    out = tmp_path / "topo.data"
    topo_data.write(out_file=out)

    assert _field_names(out) == ["topo_missing", "test_topography",
                                 "ntopofiles", "override_order",
                                 "topo_path", "topo_type",
                                 "topo_path", "topo_type"]
    # ntopofiles must still count the files, not the records.
    assert _payload(out)[2].split()[0] == "2"


@pytest.mark.python
def test_topo_data_2d_layout_unchanged_by_1d_support(tmp_path):
    r"""The default (2D) topo.data still carries the full block."""
    import clawpack.geoclaw.topotools as topotools

    topo = topotools.Topography()
    topo.path = "topo.tt3"
    topo.topo_type = 3

    topo_data = clawpack.geoclaw.data.TopographyData()
    assert topo_data.num_dim == 2, "2D must remain the default"
    topo_data.topofiles = [topo]
    out = tmp_path / "topo.data"
    topo_data.write(out_file=out)

    text = Path(out).read_text()
    assert "override_order" not in text
    for attr in ("crop_extent", "coarsen", "buffer", "align",
                 "x_shift", "y_shift", "z_shift", "negate_z"):
        assert attr in text


@pytest.mark.python
def test_dtopo_data_1d_path_line_has_no_trailing_comment(tmp_path):
    r"""1D dtopo.data path line must be bare, and carry no preprocessing block.

    ``read(iunit,*) dtopofname, dtopotype`` spans records, so the type may sit
    on the following line -- but a trailing comment on the path line is
    consumed as item 2 and aborts with "Bad integer for item 2 in list input".
    """
    import clawpack.geoclaw.dtopotools as dtopotools

    d = dtopotools.DTopography()
    d.path = "dtopo_okada.dtt1"
    d.dtopo_type = 1

    dtopo_data = clawpack.geoclaw.data.DTopoData(num_dim=1)
    dtopo_data.dtopofiles = [d]
    dtopo_data.dt_max_dtopo = 0.5
    out = tmp_path / "dtopo.data"
    dtopo_data.write(out_file=out)

    lines = _payload(out)
    # mdtopofiles, path, dtopo_type, dt_max_dtopo
    assert len(lines) == 4, lines
    assert lines[1].rstrip().endswith("'"), (
        "path line must not carry a trailing comment: %r" % lines[1])
    assert "dtopo_type" in lines[2]
    assert "dt_max_dtopo" in lines[3]
    assert "crop_extent" not in Path(out).read_text()


@pytest.mark.python
def test_1d_preprocessing_request_warns(tmp_path):
    r"""Preprocessing asked for in 1D is reported, not silently dropped."""
    import clawpack.geoclaw.topotools as topotools

    topo = topotools.Topography()
    topo.path = "celledges.data"
    topo.topo_type = 1
    topo.coarsen = 4
    topo.z_shift = 2.0

    topo_data = clawpack.geoclaw.data.TopographyData(num_dim=1)
    topo_data.topofiles = [topo]

    with pytest.warns(UserWarning, match="not supported for 1D"):
        topo_data.write(out_file=tmp_path / "topo.data")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
