#!/usr/bin/env python
# encoding: utf-8

"""Unit tests for the meteorological-forcing object model (Phase 1).

These tests exercise the new object model introduced by the met-forcing
refactor -- :class:`Track`, :class:`StormTrack`,
:class:`ParametricMetForcing`, and :class:`GriddedMetForcing` -- directly
(i.e. not only through the ``Storm`` compatibility wrapper).  They assert:

* the four objects construct;
* the readers (``read_geoclaw`` / ``read_atcf`` / ``read_data``) return the new
  objects;
* a new-object ``write_geoclaw`` / ``write_data`` produces bytes byte-identical
  to the ``Storm`` path *and* to the committed Phase-0 write goldens in
  ``tests/data/storm/characterization/``.
"""

from pathlib import Path
import sys

import numpy as np
import pytest

import clawpack.geoclaw.surge.storm as storm
from clawpack.geoclaw.surge.track import Track, StormTrack
from clawpack.geoclaw.surge.parametric import ParametricMetForcing
from clawpack.geoclaw.surge.gridded import GriddedMetForcing

# ``tests/`` has no package __init__; make the sibling helpers importable.
sys.path.insert(0, str(Path(__file__).parent))
from test_storm import _storm_input_path, _storm_check_path  # noqa: E402
from test_storm_characterization import _descriptor_head, golden_dir  # noqa: E402


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_objects_construct():
    """All four core objects construct with sensible empty defaults."""
    track = Track()
    assert track.t is None
    assert track.center is None
    # eye_location is an alias of center; ID an alias of id.
    track.eye_location = np.zeros((2, 2))
    assert track.center is track.eye_location

    storm_track = StormTrack()
    assert isinstance(storm_track, Track)
    for field in ("max_wind_speed", "max_wind_radius", "central_pressure",
                  "storm_radius", "classification", "basin", "wind_speeds"):
        assert getattr(storm_track, field) is None

    parametric = ParametricMetForcing()
    assert isinstance(parametric.track, StormTrack)
    assert parametric.file_paths == []

    gridded = GriddedMetForcing()
    assert gridded.scaling == [1.0, 1.0]
    assert gridded.crop_extent is None
    assert gridded.met_variable_map == {}


# ---------------------------------------------------------------------------
# Readers produce the new objects
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_read_geoclaw_returns_parametric():
    """ParametricMetForcing.read_geoclaw returns a ParametricMetForcing view."""
    forcing = ParametricMetForcing.read_geoclaw(_storm_check_path("tcvitals"))
    assert isinstance(forcing, ParametricMetForcing)
    assert isinstance(forcing.track, StormTrack)
    assert forcing.file_format == "geoclaw"
    assert forcing.eye_location.shape[1] == 2
    assert forcing.t.dtype.kind == "M"  # datetime64 track axis


@pytest.mark.python
@pytest.mark.storm
def test_read_atcf_returns_stormtrack():
    """StormTrack.read_atcf returns a populated StormTrack."""
    pytest.importorskip("pandas")
    track = StormTrack.read_atcf(_storm_input_path("atcf"))
    assert isinstance(track, StormTrack)
    assert track.file_format == "atcf"
    assert track.basin == "Atlantic"
    assert track.center.shape[1] == 2
    assert track.wind_speeds is not None


@pytest.mark.python
@pytest.mark.storm
def test_read_data_returns_gridded(tmp_path):
    """GriddedMetForcing.read_data returns a populated GriddedMetForcing."""
    # Produce an ASCII descriptor with the new object, then read it back.
    gridded = _make_ascii_gridded()
    descriptor = tmp_path / "ascii.storm"
    gridded.write_data(descriptor)

    read = GriddedMetForcing.read_data(descriptor)
    assert isinstance(read, GriddedMetForcing)
    assert read.file_format == 1
    assert read.crop_extent == [-100.0, -60.0, 10.0, 40.0]
    assert read.ramp_width == 3.0
    assert read.x_shift == 1.25
    assert read.y_shift == -0.5
    assert len(read.file_paths) == 4


# ---------------------------------------------------------------------------
# Write byte-equality: new object == Storm path == committed golden
# ---------------------------------------------------------------------------

WRITE_GEOCLAW_FORMATS = ["atcf", "tcvitals"]


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format", WRITE_GEOCLAW_FORMATS)
def test_write_geoclaw_matches_storm_and_golden(tmp_path, file_format):
    """New-object write_geoclaw == Storm-path write == committed golden."""
    if file_format == "atcf":
        pytest.importorskip("pandas")

    input_path = _storm_input_path(file_format)
    reader = getattr(StormTrack, "read_%s" % file_format)

    # New-object path.
    track = reader(input_path)
    forcing = ParametricMetForcing(track=track)
    new_out = tmp_path / f"{file_format}_new.storm"
    forcing.write_geoclaw(new_out)

    # Storm compatibility path.
    s = storm.Storm(input_path, file_format=file_format)
    storm_out = tmp_path / f"{file_format}_storm.storm"
    s.write(storm_out, file_format="geoclaw")

    new_bytes = new_out.read_text()
    golden = (golden_dir / f"write_geoclaw_{file_format}.txt").read_text()

    assert new_bytes == storm_out.read_text()
    assert new_bytes == golden


def _make_ascii_gridded():
    """A deterministic OWI/NWS12 ASCII GriddedMetForcing (fixed controls)."""
    gridded = GriddedMetForcing()
    gridded.time_offset = np.datetime64("2012-08-29")
    gridded.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    gridded.ramp_width = 3
    gridded.x_shift = 1.25
    gridded.y_shift = -0.5
    gridded.file_format = "ascii"
    gridded.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                          Path("storm_2.PRE"), Path("storm_2.WIN")]
    return gridded


@pytest.mark.python
@pytest.mark.storm
def test_write_data_matches_storm_and_golden(tmp_path):
    """New-object write_data == Storm-path write == committed golden (ASCII)."""
    # New-object path.
    new_desc = tmp_path / "ascii_new.storm"
    _make_ascii_gridded().write_data(new_desc)

    # Storm compatibility path.
    s = storm.Storm()
    s.time_offset = np.datetime64("2012-08-29")
    s.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    s.ramp_width = 3
    s.x_shift = 1.25
    s.y_shift = -0.5
    s.file_format = "ascii"
    s.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                    Path("storm_2.PRE"), Path("storm_2.WIN")]
    storm_desc = tmp_path / "ascii_storm.storm"
    s.write(storm_desc, file_format="data")

    new_head = _descriptor_head(new_desc.read_text())
    golden = (golden_dir / "write_data_ascii.txt").read_text()

    assert new_head == _descriptor_head(storm_desc.read_text())
    assert new_head == golden


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.netcdf
def test_write_data_netcdf_matches_golden(tmp_path):
    """New-object write_data (NWS13 netCDF) == committed descriptor golden."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    from test_storm import create_nws13_storm_file

    nc = tmp_path / "nws13.nc"
    create_nws13_storm_file(nc)

    gridded = GriddedMetForcing()
    gridded.time_offset = np.datetime64("2012-08-29")
    gridded.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    gridded.ramp_width = 3
    gridded.x_shift = 1.25
    gridded.y_shift = -0.5
    gridded.file_format = "nws13"
    gridded.file_paths = [nc]
    descriptor = tmp_path / "nws13.storm"
    gridded.write_data(descriptor,
                       var_mapping={"wind_u": "uwnd", "wind_v": "vwnd",
                                    "pressure": "press"})

    new_head = _descriptor_head(descriptor.read_text())
    golden = (golden_dir / "write_data_nws13.txt").read_text()
    assert new_head == golden


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
