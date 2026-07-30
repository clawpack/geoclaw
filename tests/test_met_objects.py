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

import clawpack.geoclaw.met.storm as storm
from clawpack.geoclaw.met.track import Track, StormTrack
from clawpack.geoclaw.met.parametric import ParametricMetForcing
from clawpack.geoclaw.met.gridded import GriddedMetForcing

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


# ---------------------------------------------------------------------------
# Regression coverage for previously-latent surge helper bugs (wrap-up fixes)
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_storm_str_datetime64():
    """Storm.__str__ renders a datetime64 track axis without error.

    Regression: the datetime64 branch used a typo (``np.datetiem64``) and
    ``.isoformat()``, which np.datetime64 does not provide.
    """
    s = storm.Storm()
    s.name = "TESTSTORM"
    s.t = np.array(["2020-08-01T00:00", "2020-08-01T06:00"],
                   dtype="datetime64[s]")
    s.file_paths = ["a.storm"]
    text = str(s)
    assert "TESTSTORM" in text
    assert "2020-08-01T00:00:00" in text
    assert "2020-08-01T06:00:00" in text


@pytest.mark.python
@pytest.mark.storm
def test_construct_fields_resolves_radius():
    """construct_fields resolves its radius argument (no NameError) and reaches
    the not-yet-implemented model stub.

    Regression: the call passed an undefined name ``x`` instead of ``r``.
    """
    s = storm.Storm()
    with pytest.raises(NotImplementedError):
        storm.construct_fields(s, 1.0e3, 0.0, model="holland_1980")


@pytest.mark.python
@pytest.mark.storm
def test_make_multi_structure_splits_by_storm(tmp_path):
    """make_multi_structure splits a multi-storm ATCF into per-storm Storms.

    Regression: the helper referenced ``os`` without importing it, and grouped
    by timestamp rather than storm identity (colliding on os.mkdir).
    """
    pytest.importorskip("pandas")
    from clawpack.geoclaw.met.tools import make_multi_structure

    # Build a 2-storm fixture from real ATCF records: relabel the cyclone
    # number (field 1) so the two synthetic storms are distinguishable.
    with open(_storm_input_path("atcf")) as data_file:
        records = [line for line in data_file
                   if len(line.split(",")) > 8][:4]
    assert records, "expected ATCF records in the test fixture"
    basin = records[0].split(",")[0].strip()

    def relabel(line, cyclone):
        fields = line.split(",")
        fields[1] = " %s" % cyclone
        return ",".join(fields)

    fixture = tmp_path / "multi.atcf"
    fixture.write_text("".join([relabel(r, "09") for r in records]
                               + [relabel(r, "11") for r in records]))

    storms = make_multi_structure(str(fixture),
                                  output_dir=str(tmp_path / "clipped"))
    assert list(storms.keys()) == [basin + "09", basin + "11"]
    for split in storms.values():
        assert isinstance(split, storm.Storm)
        assert len(split.t) == len(records)


@pytest.mark.python
@pytest.mark.storm
def test_surgedata_lives_in_geoclaw_data():
    """SurgeData/FrictionData live in clawpack.geoclaw.data, not surge.data.

    Regression: isaac/setplot_kml.py imported the nonexistent
    ``clawpack.geoclaw.surge.data`` module.
    """
    import clawpack.geoclaw.data as geodata
    assert hasattr(geodata, "SurgeData")
    assert hasattr(geodata, "FrictionData")
    with pytest.raises(ImportError):
        import clawpack.geoclaw.surge.data  # noqa: F401


@pytest.mark.python
@pytest.mark.storm
def test_isaac_setplot_kml_imports():
    """isaac/setplot_kml.py imports cleanly after the surge.data import fix."""
    pytest.importorskip("matplotlib")
    pytest.importorskip("clawpack.visclaw")
    import importlib.util

    kml_path = (Path(__file__).parents[1] / "examples" / "storm-surge"
                / "isaac" / "setplot_kml.py")
    if not kml_path.exists():
        pytest.skip("isaac/setplot_kml.py not present")
    spec = importlib.util.spec_from_file_location("isaac_setplot_kml", kml_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert hasattr(module, "setplot")


# ---------------------------------------------------------------------------
# met / surge package aliasing (Decision 1 -> B)
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_met_package_reexports_model():
    """``clawpack.geoclaw.met`` re-exports the storm/met object model."""
    import clawpack.geoclaw.met as met
    for name in ("Storm", "Track", "StormTrack", "ParametricMetForcing",
                 "GriddedMetForcing", "OWIData", "SurgeData", "MetData",
                 "construct_fields", "available_formats", "available_models"):
        assert hasattr(met, name), f"clawpack.geoclaw.met is missing {name}"


@pytest.mark.python
@pytest.mark.storm
def test_metdata_is_surgedata_alias():
    """``clawpack.geoclaw.data`` exposes ``MetData`` as an alias of ``SurgeData``."""
    import clawpack.geoclaw.data as geodata
    assert geodata.MetData is geodata.SurgeData


@pytest.mark.python
@pytest.mark.storm
def test_surge_import_emits_deprecation_warning():
    """Importing the deprecated ``surge`` package warns and points at ``met``."""
    import importlib
    import warnings
    import clawpack.geoclaw.surge as surge
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        importlib.reload(surge)   # re-run the package body to observe the warning
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert len(dep) >= 1, "surge import did not emit a DeprecationWarning"
    assert "clawpack.geoclaw.met" in str(dep[0].message)


@pytest.mark.python
@pytest.mark.storm
def test_surge_shim_delegates_to_met():
    """``surge.*`` names resolve to the identical objects in ``met.*``."""
    import clawpack.geoclaw.met.storm as met_storm
    import clawpack.geoclaw.met.track as met_track
    import clawpack.geoclaw.surge.storm as surge_storm
    import clawpack.geoclaw.surge.track as surge_track
    assert surge_storm.Storm is met_storm.Storm
    assert surge_track.StormTrack is met_track.StormTrack
    # Private names delegate too (gridded.py imports _Meta from track).
    assert surge_track._Meta is met_track._Meta


# ---------------------------------------------------------------------------
# OWI (Oceanweather WIN/PRE) field I/O
# ---------------------------------------------------------------------------

def _make_owi_data(nt=3, ny=4, nx=5):
    """Build a small deterministic OWIData object."""
    from clawpack.geoclaw.met.data_storms import OWIData

    lon = -99.0 + np.arange(nx) * 0.25
    lat = 8.0 + np.arange(ny) * 0.25
    t = (np.datetime64("2012-08-20T12:00:00")
         + (np.arange(nt) * 3600).astype("timedelta64[s]"))
    ramp = np.arange(nt * ny * nx, dtype=float).reshape(nt, ny, nx)
    return OWIData(time=t, longitude=lon, latitude=lat,
                   wind_u=ramp, wind_v=-ramp, pressure=101300.0 - ramp)


@pytest.mark.python
@pytest.mark.storm
def test_owi_roundtrip(tmp_path):
    """write_owi -> read_owi reproduces every field exactly."""
    from clawpack.geoclaw.met import data_storms

    d = _make_owi_data()
    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    data_storms.write_owi(d, pre, win)
    r = data_storms.read_owi(pre, win)

    assert np.array_equal(d.time, r.time)
    np.testing.assert_allclose(d.longitude, r.longitude)
    np.testing.assert_allclose(d.latitude, r.latitude)
    np.testing.assert_allclose(d.wind_u, r.wind_u)
    np.testing.assert_allclose(d.wind_v, r.wind_v)
    # Pressure survives the Pa<->mbar conversion to the f10.4 precision.
    np.testing.assert_allclose(d.pressure, r.pressure, atol=1e-2)
    # The cheap start-time peek agrees with the first record.
    assert data_storms.read_owi_start_time(pre) == d.time[0]


@pytest.mark.python
@pytest.mark.storm
def test_owi_written_format_is_80_col(tmp_path):
    """Written OWI lines match the fixed 80-column contract the Fortran reads."""
    from clawpack.geoclaw.met import data_storms

    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    data_storms.write_owi(_make_owi_data(), pre, win)
    lines = win.read_text().splitlines()
    assert len(lines[0]) == 80                       # file header
    assert lines[0].startswith("Oceanweather WIN/PRE Format")
    assert len(lines[1]) == 80                       # grid-spec header
    assert lines[1].startswith("iLat=") and "iLong=" in lines[1]
    assert "DT=" in lines[1]


@pytest.mark.python
@pytest.mark.storm
def test_owi_read_real_isaac():
    """The reader parses the committed real isaac.WIN/PRE sample."""
    from clawpack.geoclaw.met import data_storms

    isaac = (Path(__file__).parents[1] / "examples" / "storm-surge" / "isaac")
    if not (isaac / "isaac.PRE").exists():
        pytest.skip("isaac OWI sample not present")

    d = data_storms.read_owi(isaac / "isaac.PRE", isaac / "isaac.WIN")
    assert d.shape == (96, 116)                       # (ny, nx)
    assert d.wind_u.shape == (d.num_times, 96, 116)
    np.testing.assert_allclose(d.longitude[0], -99.0)
    np.testing.assert_allclose(d.latitude[0], 8.0)
    # Ambient far-field pressure: 1013 mbar -> 101300 Pa.
    np.testing.assert_allclose(d.pressure[0].max(), 101300.0)


@pytest.mark.python
@pytest.mark.storm
def test_gridded_from_owi_descriptor(tmp_path):
    """from_owi builds a format-1 descriptor pointing at the WIN/PRE pair."""
    from clawpack.geoclaw.met import data_storms

    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    data_storms.write_owi(_make_owi_data(), pre, win)

    forcing = GriddedMetForcing.from_owi(pre, win)
    assert forcing.file_format == "owi"
    assert forcing.file_paths == [pre, win]           # [PRE, WIN] order
    assert forcing.time_offset == np.datetime64("2012-08-20T12:00:00")

    # A format-1 (ascii/owi) descriptor round-trips through read_data.
    descriptor = tmp_path / "met.storm"
    forcing.write_data(descriptor)
    back = GriddedMetForcing.read_data(descriptor)
    assert back.file_format == 1
    assert [p.name for p in back.file_paths] == ["x.PRE", "x.WIN"]


@pytest.mark.python
@pytest.mark.storm
def test_gridded_to_owi(tmp_path):
    """to_owi writes the pair and points the forcing at it."""
    from clawpack.geoclaw.met import data_storms

    d = _make_owi_data()
    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    forcing = GriddedMetForcing().to_owi(d, pre, win)

    assert pre.exists() and win.exists()
    assert forcing.file_format == "owi"
    assert forcing.file_paths == [pre, win]
    # The files it wrote read back to the same fields.
    r = data_storms.read_owi(pre, win)
    np.testing.assert_allclose(d.wind_u, r.wind_u)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
