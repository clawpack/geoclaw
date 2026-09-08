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
import warnings

import numpy as np
import pytest

import clawpack.geoclaw.met.storm as storm
from clawpack.geoclaw.met.track import Track, StormTrack, fill_rad_w_other_source
from clawpack.geoclaw.met.parametric import ParametricMetForcing
from clawpack.geoclaw.met.gridded import GriddedMetForcing

# ``tests/`` has no package __init__; make the sibling helpers importable.
sys.path.insert(0, str(Path(__file__).parent))
from test_storm import (_storm_input_path, _storm_check_path,  # noqa: E402
                        _make_storm_from_format)
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


# ---------------------------------------------------------------------------
# IBTrACS reader output types.
#
# read_ibtracs is the only reader backed by xarray, and it used to let xarray
# types leak out: ``t`` stayed a DataArray and ``classification`` stayed bytes.
# These are contract tests rather than snapshots, so they keep holding as the
# fixture or the xarray version changes.
# ---------------------------------------------------------------------------

_IBTRACS_KWARGS = {"sid": "2008245N17323", "agency_pref": ["wmo", "usa"]}


@pytest.mark.python
@pytest.mark.storm
def test_ibtracs_reader_types():
    """read_ibtracs returns plain NumPy/Python types, like every other reader."""
    pytest.importorskip("xarray")
    track = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                    **_IBTRACS_KWARGS)

    # Times: a datetime64 ndarray, not an xarray DataArray.  Indexing a
    # DataArray yields 0-d DataArrays, which is what broke
    # fill_rad_w_other_source.
    assert isinstance(track.t, np.ndarray)
    assert track.t.dtype == np.dtype("datetime64[s]")
    assert isinstance(track.t[0], np.datetime64)

    # time_offset is drawn from t, so it follows.
    assert isinstance(track.time_offset, np.datetime64)

    # Classification decoded from the netCDF's bytes, so it compares against
    # ordinary string literals rather than b'TD'.
    assert track.classification.dtype.kind == "U"
    assert "TD" in track.classification
    assert not any(str(value).startswith("b'") for value in track.classification)

    assert track.event.dtype.kind == "U"


@pytest.mark.python
@pytest.mark.storm
def test_ibtracs_second_resolution_times():
    """IBTrACS times land on whole seconds, not a nanosecond roundoff tail.

    IBTrACS stores times that decode to values like
    ``2008-09-01T06:00:00.000039936``.  Carried through to ``write_geoclaw``
    those turn nominal hour offsets into ``-3.60000003e+03`` instead of
    ``-3.60000000e+03``; the committed ``ibtracs_geoclaw.txt`` baseline has the
    exact values, so truncating to seconds restores agreement with it.
    """
    pytest.importorskip("xarray")
    track = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                    **_IBTRACS_KWARGS)

    offsets = (track.t - track.time_offset) / np.timedelta64(1, "s")
    # Every IBTrACS observation is on a whole minute; no sub-second residue.
    assert np.allclose(offsets % 60.0, 0.0, atol=0.0)


@pytest.mark.python
@pytest.mark.storm
def test_ibtracs_reader_no_future_warnings():
    """read_ibtracs uses no deprecated xarray API.

    ``Dataset.dims`` returning a mapping is deprecated and becomes a set of
    dimension names; ``Dataset.sizes`` is the replacement.  Promoting the
    warning to an error keeps this from silently rotting on the next xarray
    release.
    """
    pytest.importorskip("xarray")
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        # The reader legitimately warns about missing RMW/ROCI; that is not the
        # class of warning under test here.
        warnings.filterwarnings("ignore", category=UserWarning)
        StormTrack.read_ibtracs(_storm_input_path("ibtracs"), **_IBTRACS_KWARGS)


@pytest.mark.python
@pytest.mark.storm
def test_fill_rad_accepts_scalar_or_zero_d_time():
    """fill_rad_w_other_source accepts a scalar and a 0-d DataArray time alike.

    Callers pass ``storm.t[n]``, so whichever of those a reader produces has to
    work.  A 0-d DataArray used to raise
    ``ValueError: Could not convert object to NumPy datetime``.
    """
    xr = pytest.importorskip("xarray")
    pytest.importorskip("pandas")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        target = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                         **_IBTRACS_KWARGS)
        source = StormTrack.read_atcf(_storm_input_path("atcf"))

    scalar_t = target.t[len(target.t) // 2]
    from_scalar = fill_rad_w_other_source(scalar_t, target, source,
                                          "max_wind_radius")

    zero_d = xr.DataArray(scalar_t)
    assert zero_d.ndim == 0
    from_zero_d = fill_rad_w_other_source(zero_d, target, source,
                                          "max_wind_radius")

    assert np.isclose(from_scalar, from_zero_d)


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
# Missing-data contract
#
# np.nan is the in-memory marker in every reader; write_geoclaw additionally
# tolerates a negative for callers still on the v5.9.0 -1 convention; zero is
# never missing.  These lock all three halves of that contract.
# ---------------------------------------------------------------------------

def _synthetic_track(**overrides):
    """A minimal four-point ParametricMetForcing for writer tests."""
    n = 4
    forcing = ParametricMetForcing()
    forcing.t = np.array(["2020-01-01T00", "2020-01-01T06",
                          "2020-01-01T12", "2020-01-01T18"],
                         dtype="datetime64[s]")
    forcing.time_offset = forcing.t[0]
    forcing.eye_location = np.column_stack([np.linspace(-80.0, -83.0, n),
                                            np.linspace(25.0, 28.0, n)])
    forcing.max_wind_speed = np.full(n, 40.0)
    forcing.central_pressure = np.full(n, 98000.0)
    forcing.max_wind_radius = np.full(n, 40e3)
    forcing.storm_radius = np.full(n, 300e3)
    for name, value in overrides.items():
        setattr(forcing, name, value)
    return forcing


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("field", ["max_wind_speed", "central_pressure",
                                   "max_wind_radius", "storm_radius"])
def test_missing_marker_equivalence(tmp_path, field):
    """A NaN-marked and a -1-marked storm write byte-identical files.

    The readers now emit NaN, but objects built by hand and community fill
    functions may still use -1.  Both must reach the same fill/skip branch.
    """
    nan_values = np.full(4, np.nan)
    sentinel_values = np.full(4, -1.0)

    nan_path = tmp_path / "nan.storm"
    sentinel_path = tmp_path / "sentinel.storm"

    _synthetic_track(**{field: nan_values}).write_geoclaw(nan_path)
    with warnings.catch_warnings():
        # The -1 path warns that the convention is deprecated; that is the
        # documented behavior, not a failure.
        warnings.simplefilter("ignore", DeprecationWarning)
        _synthetic_track(**{field: sentinel_values}).write_geoclaw(
            sentinel_path)

    assert nan_path.read_bytes() == sentinel_path.read_bytes()


@pytest.mark.python
@pytest.mark.storm
def test_failed_fill_is_skipped_not_written(tmp_path):
    """A fill that cannot produce a value leaves the row skipped.

    Returning NaN or -1 from a fill used to be written into the file verbatim,
    producing a storm file GeoClaw cannot run.
    """
    for sentinel in (np.nan, -1.0):
        path = tmp_path / f"fill_{sentinel}.storm"
        forcing = _synthetic_track(max_wind_radius=np.full(4, np.nan))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            forcing.write_geoclaw(
                path, fill_dict={"max_wind_radius": lambda t, s: sentinel})
        # Header line is the cast count; every row should have been skipped.
        assert path.read_text().splitlines()[0].strip() == "0"


@pytest.mark.python
@pytest.mark.storm
def test_zero_is_not_missing(tmp_path):
    """Zero is real data, not a missing marker.

    HURDAT2 reports genuinely-zero quadrant wind radii for weak systems.  A
    'treat <= 0 as missing' rule would silently reinterpret those.
    """
    path = tmp_path / "zero.storm"
    forcing = _synthetic_track(max_wind_speed=np.zeros(4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # No fill supplied: if zero were treated as missing every row would be
        # skipped and this would be a zero-cast file.
        forcing.write_geoclaw(path)

    lines = path.read_text().splitlines()
    assert lines[0].strip() == "4"
    assert float(lines[3].split()[3]) == 0.0


@pytest.mark.python
@pytest.mark.storm
def test_write_geoclaw_does_not_mutate_storm(tmp_path):
    """Writing twice with different fills honors the second fill.

    Fills used to be written back onto the storm object, so the second write saw
    a storm that was no longer missing and silently reused the first fill.
    """
    forcing = _synthetic_track(max_wind_radius=np.full(4, np.nan))

    first = tmp_path / "first.storm"
    second = tmp_path / "second.storm"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        forcing.write_geoclaw(
            first, fill_dict={"max_wind_radius": lambda t, s: 20e3})
        forcing.write_geoclaw(
            second, fill_dict={"max_wind_radius": lambda t, s: 60e3})

    assert float(first.read_text().splitlines()[3].split()[4]) == 20e3
    assert float(second.read_text().splitlines()[3].split()[4]) == 60e3
    # And the storm itself is untouched.
    assert np.all(np.isnan(forcing.max_wind_radius))


@pytest.mark.python
@pytest.mark.storm
def test_fill_dict_storm_radius_not_clobbered(tmp_path):
    """A caller-supplied storm_radius fill beats the built-in 500 km default.

    The default was applied with dict.update() *after* copying the caller's
    fill_dict, so a caller-supplied storm_radius fill was silently discarded.
    """
    path = tmp_path / "roci.storm"
    forcing = _synthetic_track(storm_radius=np.full(4, np.nan))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        forcing.write_geoclaw(
            path, fill_dict={"storm_radius": lambda t, s: 250e3})

    assert float(path.read_text().splitlines()[3].split()[6]) == 250e3


@pytest.mark.python
@pytest.mark.storm
def test_hurdat_sentinels_become_nan(tmp_path):
    """HURDAT2 -99 wind / -999 pressure become NaN, not scaled sentinels.

    Converting units before normalizing turned -999 mbar into -99900.0 Pa, a
    value that looks physical enough to be written out.
    """
    header = "AL011980,          UNNAMED,      2,\n"
    rows = ("19800101, 0000,  , TD, 25.0N,  80.0W, -99, -999,"
            + " -999," * 11 + " -999,\n"
            "19800101, 0600,  , TD, 25.5N,  80.5W,  35, 1005,"
            + " -999," * 11 + " -999,\n")
    path = tmp_path / "sentinel_hurdat.txt"
    path.write_text(header + rows)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        track = StormTrack.read_hurdat(path)

    assert np.isnan(track.max_wind_speed[0])
    assert np.isnan(track.central_pressure[0])
    # The valid second row is untouched and correctly converted.
    assert np.isclose(track.central_pressure[1], 100500.0)
    assert not np.isnan(track.max_wind_speed[1])


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format", ["atcf", "tcvitals", "ibtracs"])
def test_written_radii_are_positive(tmp_path, file_format):
    """Every storm file written from a bundled input is runnable.

    A non-positive storm_radius zeros the forcing through the spatial ramp and a
    non-positive max_wind_radius divides by zero in the Holland profiles, so a
    file with either is not a valid GeoClaw input.  hurdat and jma are excluded
    until met.reconstruction supplies real RMW/ROCI estimators; their current
    baselines use placeholder zeros.
    """
    if file_format == "ibtracs":
        pytest.importorskip("xarray")
    pytest.importorskip("pandas")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
        write_kwargs = {}
        if fill_mwr is not None:
            write_kwargs["max_wind_radius_fill"] = fill_mwr
        if fill_rad is not None:
            write_kwargs["storm_radius_fill"] = fill_rad
        path = tmp_path / f"{file_format}.storm"
        test_storm.write(path, file_format="geoclaw", **write_kwargs)

    values = np.loadtxt(path, skiprows=3)
    assert values.shape[0] > 0
    assert np.all(np.isfinite(values))
    assert np.all(values[:, 4] > 0.0), "max_wind_radius must be positive"
    assert np.all(values[:, 6] > 0.0), "storm_radius must be positive"


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


# ---------------------------------------------------------------------------
# Multi-storm archive ingestion
#
# Both HURDAT2 and IBTrACS ship one file per basin holding every storm on
# record.  The `iter_*` generators and the `read_*` classmethods share a single
# parser, so the load-bearing check is that the two paths produce identical
# objects -- not two independent goldens that could drift apart.
# ---------------------------------------------------------------------------

def _parsed_state(track):
    """Canonical text of a track's state, excluding provenance.

    Compared as serialized text rather than as a dict, because the state
    contains NaN and ``nan != nan`` makes dict equality fail spuriously -- the
    same reason the characterization tests compare golden *text*.
    ``file_paths`` is dropped: it records where the object came from, which
    legitimately differs between a bulk walk and a single-storm read, while
    everything else is the parsed track and must match exactly.
    """
    from test_storm_characterization import _storm_state, _dumps
    state = _storm_state(track)
    state.pop("file_paths", None)
    return _dumps(state)


def _two_storm_hurdat(tmp_path):
    """A two-storm HURDAT2 file, relabeled from the bundled single-storm one."""
    original = _storm_input_path("hurdat").read_text().rstrip("\n").split("\n")
    header, records = original[0], original[1:]

    # Re-frame with an accurate declared record count, then a second storm with
    # a different id and name.
    first = f"AL082008,              IKE,{len(records):>4},"
    second = f"AL102008,             JOSE,{len(records):>4},"
    path = tmp_path / "two_storms.txt"
    path.write_text("\n".join([first] + records + [second] + records) + "\n")
    assert header  # the bundled fixture really does start with a header
    return path


@pytest.mark.python
@pytest.mark.storm
def test_hurdat_declared_record_count_is_enforced(tmp_path):
    """A HURDAT2 header declares its record count; a mismatch must be caught.

    This is the framing invariant the format hands us for free.  Without it the
    single-storm reader walked into the *next* storm's header and fed
    ``"AL092008"`` to ``np.datetime64``.
    """
    path = _two_storm_hurdat(tmp_path)
    lines = path.read_text().split("\n")

    # Truncate the first storm's declared count by one.
    from clawpack.geoclaw.met.track import _hurdat_blocks
    blocks = list(_hurdat_blocks(path))
    assert len(blocks) == 2
    assert all(len(records) == header["num_records"]
               for header, records in blocks)
    total = sum(header["num_records"] for header, _ in blocks)
    data_lines = sum(1 for line in lines
                     if line[:8].isdigit())
    assert total == data_lines, "declared counts must sum to the data lines"

    # A header claiming more records than follow is an error, not a silent read
    # into the next storm.
    lines[0] = lines[0].replace(f"{len(blocks[0][1]):>4},", "  99,")
    bad = tmp_path / "bad_count.txt"
    bad.write_text("\n".join(lines))
    with pytest.raises(ValueError, match="declares"):
        list(_hurdat_blocks(bad))


@pytest.mark.python
@pytest.mark.storm
def test_iter_hurdat_matches_read_hurdat(tmp_path):
    """Bulk and single-storm HURDAT2 reads produce identical objects."""
    from clawpack.geoclaw.met.track import iter_hurdat

    path = _two_storm_hurdat(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tracks = list(iter_hurdat(path))
        selected = StormTrack.read_hurdat(path, storm_id="AL102008")

    assert [t.ID for t in tracks] == ["AL082008", "AL102008"]
    assert [t.name for t in tracks] == ["IKE", "JOSE"]
    assert _parsed_state(tracks[1]) == _parsed_state(selected)


@pytest.mark.python
@pytest.mark.storm
def test_read_hurdat_requires_a_selector_for_multi_storm_files(tmp_path):
    """Ambiguity is an error naming the candidates, not an arbitrary pick."""
    path = _two_storm_hurdat(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="2 storms"):
            StormTrack.read_hurdat(path)
        with pytest.raises(ValueError, match="No storm"):
            StormTrack.read_hurdat(path, name="NOSUCHSTORM")
        # A single-storm file still needs no selector.
        single = StormTrack.read_hurdat(_storm_input_path("hurdat"))
        assert single.ID == "AL082008"


@pytest.mark.python
@pytest.mark.storm
def test_iter_hurdat_filters(tmp_path):
    """Filters combine, and years accepts an int or an iterable."""
    from clawpack.geoclaw.met.track import iter_hurdat, catalog_hurdat

    path = _two_storm_hurdat(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert len(list(iter_hurdat(path, names=["JOSE"]))) == 1
        assert len(list(iter_hurdat(path, years=2008))) == 2
        assert len(list(iter_hurdat(path, years=range(2000, 2005)))) == 0
        assert len(list(iter_hurdat(path, basins="AL"))) == 2
        assert len(list(iter_hurdat(path, basins="EP"))) == 0

    catalog = catalog_hurdat(path)
    assert list(catalog.storm_id) == ["AL082008", "AL102008"]
    assert set(catalog.basin) == {"Atlantic"}
    assert (catalog.t_start <= catalog.t_end).all()


@pytest.mark.python
@pytest.mark.storm
def test_iter_ibtracs_matches_read_ibtracs(tmp_path):
    """Bulk and single-storm IBTrACS reads produce identical objects."""
    xr = pytest.importorskip("xarray")
    from clawpack.geoclaw.met.track import iter_ibtracs

    path = _storm_input_path("ibtracs")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bulk = list(iter_ibtracs(path, agency_pref=["wmo", "usa"]))
        single = StormTrack.read_ibtracs(path, **_IBTRACS_KWARGS)

    assert len(bulk) == 1
    assert _parsed_state(bulk[0]) == _parsed_state(single)


@pytest.mark.python
@pytest.mark.storm
def test_iter_ibtracs_over_multiple_storms(tmp_path):
    """Two storms in one file are read independently, with no cross-contamination.

    The bundled fixture is a stripped single-storm slice, so a second storm is
    synthesized by relabeling a copy.  Note this cannot exercise real
    heterogeneity (differing agencies, basins, valid ranges) -- only the
    remote-marked test over a real basin file does that.
    """
    xr = pytest.importorskip("xarray")
    from clawpack.geoclaw.met.track import iter_ibtracs, catalog_ibtracs

    with xr.open_dataset(_storm_input_path("ibtracs")) as ds:
        second = ds.copy(deep=True)
        second["sid"] = second.sid.copy(data=[b"2008245N17999"])
        second["name"] = second.name.copy(data=[b"NOTIKE"])
        combined = xr.concat([ds, second], dim="storm")
        path = tmp_path / "two_storms.nc"
        combined.to_netcdf(path)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tracks = list(iter_ibtracs(path, agency_pref=["wmo", "usa"]))

    assert [t.ID for t in tracks] == ["2008245N17323", "2008245N17999"]
    assert [t.name for t in tracks] == ["IKE", "NOTIKE"]
    # Each track's times stay inside its own record -- no bleed between storms.
    for track in tracks:
        assert len(track.t) == len(tracks[0].t)
        assert np.array_equal(track.t, tracks[0].t)

    # Selecting by sid picks exactly one of them.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        agencies = {"agency_pref": ["wmo", "usa"]}
        assert len(list(iter_ibtracs(path, sids=["2008245N17999"],
                                     **agencies))) == 1
        assert len(list(iter_ibtracs(path, names=["IKE"], **agencies))) == 1

    catalog = catalog_ibtracs(path)
    assert list(catalog.sid) == ["2008245N17323", "2008245N17999"]
    assert (catalog.num_records > 0).all()


@pytest.mark.python
@pytest.mark.storm
def test_iter_ibtracs_skip_invalid(tmp_path):
    """A storm with no coincident wind and pressure is skipped, not fatal.

    Over a whole basin these exist (57 of 741 for the North Atlantic
    1980-2025), so one of them must not abort the sweep.
    """
    xr = pytest.importorskip("xarray")
    from clawpack.geoclaw.met.track import iter_ibtracs

    with xr.open_dataset(_storm_input_path("ibtracs")) as ds:
        broken = ds.copy(deep=True)
        broken["sid"] = broken.sid.copy(data=[b"2008245N17999"])
        for name in list(broken.data_vars):
            if name.endswith("_wind"):
                broken[name] = broken[name] * np.nan
        combined = xr.concat([ds, broken], dim="storm")
        path = tmp_path / "one_broken.nc"
        combined.to_netcdf(path)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tracks = list(iter_ibtracs(path, agency_pref=["wmo", "usa"]))
    assert [t.ID for t in tracks] == ["2008245N17323"]
    assert any("Skipped 1" in str(w.message) for w in caught)

    with pytest.raises((ValueError, RuntimeError)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            list(iter_ibtracs(path, agency_pref=["wmo", "usa"],
                              skip_invalid=False))


@pytest.mark.python
@pytest.mark.storm
def test_iter_atcf_matches_make_multi_structure(tmp_path):
    """The in-memory ATCF walk agrees with the split-to-disk one.

    ``make_multi_structure`` is now a thin wrapper over ``iter_atcf``; this
    pins that they stay equivalent, and that ``iter_atcf`` leaves nothing behind.
    """
    pytest.importorskip("pandas")
    from clawpack.geoclaw.met.tools import iter_atcf, make_multi_structure

    path = _storm_input_path("atcf")
    before = set(Path.cwd().iterdir())
    walked = list(iter_atcf(path))
    assert set(Path.cwd().iterdir()) == before, "iter_atcf must not write files"

    split = make_multi_structure(path, output_dir=str(tmp_path / "clipped"))
    assert [storm_id for storm_id, _ in walked] == list(split.keys())
    for storm_id, storm in walked:
        assert _parsed_state(storm) == _parsed_state(split[storm_id])


@pytest.mark.remote
@pytest.mark.storm
@pytest.mark.slow
def test_full_basin_sweep():
    """End-to-end over the real archives; the only test of real heterogeneity.

    The bundled fixtures are a single-storm HURDAT2 excerpt and a stripped
    single-storm IBTrACS slice, so nothing offline exercises differing agencies,
    basin crossings, or the ~8% of storms with no coincident wind and pressure.
    Requires the archives on disk; set ``GEOCLAW_TRACK_ARCHIVES`` to the
    directory holding ``IBTrACS.NA.v04r01.nc`` and a HURDAT2 Atlantic text file.
    """
    import os
    from clawpack.geoclaw.met.track import (iter_hurdat, iter_ibtracs,
                                            catalog_hurdat, catalog_ibtracs,
                                            _hurdat_blocks)

    archives = os.environ.get("GEOCLAW_TRACK_ARCHIVES")
    if archives is None:
        pytest.skip("set GEOCLAW_TRACK_ARCHIVES to the archive directory")
    hurdat_path = Path(archives) / "hurdat2_atlantic.txt"
    ibtracs_path = Path(archives) / "IBTrACS.NA.v04r01.nc"
    if not hurdat_path.exists() or not ibtracs_path.exists():
        pytest.skip(f"archives not found under {archives}")

    # --- HURDAT2: record-count conservation over the whole file -------------
    blocks = list(_hurdat_blocks(hurdat_path))
    declared = sum(header["num_records"] for header, _ in blocks)
    parsed = sum(len(records) for _, records in blocks)
    assert declared == parsed, "declared record counts must match what parsed"

    with open(hurdat_path) as handle:
        data_lines = sum(1 for line in handle if line[:8].isdigit())
    assert declared == data_lines, "and must account for every data line"

    catalog = catalog_hurdat(hurdat_path)
    assert len(catalog) == len(blocks)
    assert catalog.num_records.sum() == declared

    # --- IBTrACS: catalog agrees with the sweep, no cross-contamination -----
    years = range(1980, 2026)
    index = catalog_ibtracs(ibtracs_path, basin="NA", years=years)
    assert len(index) > 500, "the North Atlantic 1980-2025 has ~740 storms"
    spans = {row.sid: (row.t_start, row.t_end) for row in index.itertuples()}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tracks = list(iter_ibtracs(ibtracs_path, basin="NA", years=years))

    assert 0 < len(tracks) <= len(index)
    for track in tracks:
        assert track.t.dtype == np.dtype("datetime64[s]")
        # Times are monotone within a storm...
        assert (np.diff(track.t) >= np.timedelta64(0, "s")).all()
        # ...and stay inside that storm's own span from the catalog, which is
        # the check that would catch one storm's records bleeding into another.
        t_start, t_end = spans[track.ID]
        assert track.t[0] >= t_start
        assert track.t[-1] <= t_end
        assert track.basin == "NA"

    # --- Cross-archive: the same storm agrees between the two archives ------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ike_hurdat = next(iter_hurdat(hurdat_path, names=["IKE"], years=2008))
        ike_ibtracs = next(iter_ibtracs(ibtracs_path, names=["IKE"],
                                        years=2008))

    # Both archives are NHC best track for the Atlantic, so coincident times
    # must carry the same eye position to best-track precision (0.1 deg).
    shared, h_index, i_index = np.intersect1d(
        ike_hurdat.t, ike_ibtracs.t, return_indices=True)
    assert len(shared) > 20, "the two archives should share most of Ike's track"
    assert np.allclose(ike_hurdat.eye_location[h_index],
                       ike_ibtracs.eye_location[i_index], atol=0.11)

    # ...and on the quadrant wind radii, where both reported one.  Three
    # archives (with the ATCF fixture, offline) agreeing on all twelve radii is
    # what makes `wind_radii` trustworthy rather than merely well-shaped.
    h_radii = ike_hurdat.wind_radii[h_index]
    i_radii = ike_ibtracs.wind_radii[i_index]
    both = np.isfinite(h_radii) & np.isfinite(i_radii) & (h_radii > 0)
    assert both.sum() > 50, "expected many coincident reported radii"
    assert np.allclose(h_radii[both], i_radii[both], rtol=0, atol=1.0)

    # HURDAT2 gained a radius of maximum wind with the 2021 season.  Confirm it
    # against IBTrACS `usa_rmw`, which is the independent evidence for its
    # meaning and units (nautical miles) -- the format PDF is not machine
    # readable, so this cross-check is the provenance.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        modern = {t.name: t for t in iter_hurdat(hurdat_path,
                                                 years=range(2021, 2026))}
        compared = 0
        for track in iter_ibtracs(ibtracs_path, basin="NA",
                                  years=range(2021, 2026)):
            other = modern.get(track.name)
            if other is None:
                continue
            shared, i_ix, h_ix = np.intersect1d(track.t, other.t,
                                                return_indices=True)
            mine = track.max_wind_radius[i_ix]
            theirs = other.max_wind_radius[h_ix]
            usable = np.isfinite(mine) & np.isfinite(theirs)
            assert np.allclose(mine[usable], theirs[usable], rtol=0, atol=1.0)
            compared += int(usable.sum())
    assert compared > 1000, f"only {compared} records compared"


# ---------------------------------------------------------------------------
# Quadrant wind radii, standardized vocabularies, dtype widths
# ---------------------------------------------------------------------------

_MATURE_HURDAT = "hurdat_ike_mature.txt"


@pytest.mark.python
@pytest.mark.storm
def test_quadrant_radii_agree_across_archives():
    """ATCF and HURDAT2 report the same quadrant radii for the same storm.

    Both are NHC best track for Hurricane Ike (2008), read by two independent
    parsers, so coincident times must agree exactly.  This is the check that
    gives the new `wind_radii` field its meaning -- a shape/plumbing test would
    pass just as well with the quadrants transposed or the thresholds swapped.
    """
    pytest.importorskip("pandas")
    from clawpack.geoclaw.met.track import iter_hurdat

    atcf = StormTrack.read_atcf(_storm_input_path("atcf"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hurdat = next(iter_hurdat(_storm_input_path("hurdat").parent
                                  / _MATURE_HURDAT))

    shared, atcf_index, hurdat_index = np.intersect1d(
        atcf.t.to_numpy().astype("datetime64[s]"),
        hurdat.t, return_indices=True)
    assert len(shared) >= 8, "the fixtures should overlap in Ike's mature stage"

    a = atcf.wind_radii[atcf_index]
    h = hurdat.wind_radii[hurdat_index]
    # Compare where both reported something; ATCF maps an unavailable quadrant
    # to NaN while HURDAT2 keeps a real zero, so only the reported overlap is
    # comparable (see the StormTrack docstring).
    both = np.isfinite(a) & np.isfinite(h)
    assert both.sum() > 50, "expected many coincident reported radii"
    assert np.allclose(a[both], h[both], rtol=0, atol=1.0)


@pytest.mark.python
@pytest.mark.storm
def test_quadrant_radii_shape_and_ordering():
    """Radii shrink as the wind threshold rises, and RMW is inside r34."""
    pytest.importorskip("pandas")
    from clawpack.geoclaw.met.track import iter_hurdat

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        track = next(iter_hurdat(_storm_input_path("hurdat").parent
                                 / _MATURE_HURDAT))

    assert track.wind_radii.shape == (len(track.t), 3, 4)
    assert track.wind_radii_thresholds.shape == (3,)
    # 34 / 50 / 64 kt in m/s, increasing.
    assert (np.diff(track.wind_radii_thresholds) > 0).all()

    # A higher wind threshold cannot extend further than a lower one.
    for j in range(track.wind_radii.shape[1] - 1):
        lower, higher = track.wind_radii[:, j, :], track.wind_radii[:, j + 1, :]
        both = np.isfinite(lower) & np.isfinite(higher) & (higher > 0)
        assert (higher[both] <= lower[both] + 1.0).all()

    # Where an RMW is reported it must sit inside the 34-kt envelope.
    r34 = np.nanmean(np.where(track.wind_radii[:, 0, :] > 0,
                              track.wind_radii[:, 0, :], np.nan), axis=1)
    known = np.isfinite(track.max_wind_radius) & np.isfinite(r34)
    if known.any():
        assert (track.max_wind_radius[known] <= r34[known]).all()


@pytest.mark.python
@pytest.mark.storm
def test_hurdat_zero_radius_is_data_but_sentinel_is_not():
    """HURDAT2 distinguishes a zero radius from an unknown one; keep both.

    A weak system has no 34-kt winds, so their radius is genuinely zero -- which
    is different from ``-999``.  Locking this against a future refactor that
    treats `<= 0` as missing.
    """
    from clawpack.geoclaw.met.track import iter_hurdat

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        genesis = StormTrack.read_hurdat(_storm_input_path("hurdat"))
        mature = next(iter_hurdat(_storm_input_path("hurdat").parent
                                  / _MATURE_HURDAT))

    # The bundled fixture is Ike's 30-kt genesis: all twelve fields are 0.
    assert np.isfinite(genesis.wind_radii).all()
    assert (genesis.wind_radii == 0.0).all()

    # The mature fixture has real radii and -999 for its (pre-2021) RMW.
    assert (mature.wind_radii[:, 0, :] > 0).any()
    assert np.isnan(mature.max_wind_radius).all()


@pytest.mark.python
@pytest.mark.storm
def test_hurdat_radius_of_maximum_wind_field():
    """HURDAT2's field 21 is an RMW in nautical miles, added for 2021 onward.

    Verified against IBTrACS ``usa_rmw``: 2676 of 2676 coincident 2021-2025
    North Atlantic records agree exactly when read as nautical miles.  Older
    revisions end the line at 20 fields, which must read as missing rather than
    raising.
    """
    from clawpack.geoclaw.met.track import iter_hurdat

    # A synthetic 2021-style record with field 21 = 15 nmi.
    modern = ("AL012021,          SYNTHTC,    1,\n"
              "20210801, 0000,  , HU,  25.0N,  80.0W,  90,  960,"
              "  120,  100,   75,  120,   60,   50,   40,   50,"
              "   30,   25,   20,   25,   15,\n")
    import tempfile
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "modern.txt"
        path.write_text(modern)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            track = next(iter_hurdat(path))
    expected = 15 * 1852.0
    assert np.isclose(track.max_wind_radius[0], expected, atol=1.0)

    # The bundled fixture predates the column and must still read.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy = StormTrack.read_hurdat(_storm_input_path("hurdat"))
    assert np.isnan(legacy.max_wind_radius).all()


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format,expected", [
    ("hurdat", {"LO"}), ("atcf", {"TD", "TS", "HU", "EX"}),
    ("jma", None), ("tcvitals", None)])
def test_classification_is_not_truncated(file_format, expected):
    """``np.empty(n, dtype=str)`` is ``<U1``; multi-character codes were cut.

    HURDAT2's ``"LO"`` became ``"L"`` -- colliding with the landfall *event*
    code -- and JMA/tcvitals lost their codes the same way.  Asserted on content
    rather than dtype, since ATCF's classification comes back from pandas as
    object dtype (and, separately, with the source's leading whitespace).
    """
    pytest.importorskip("pandas")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reader = getattr(StormTrack, f"read_{file_format}")
        track = reader(_storm_input_path(file_format))

    assert track.classification is not None
    codes = {str(code).strip()
             for code in np.atleast_1d(track.classification)}
    if expected is not None:
        assert expected <= codes, f"got {sorted(codes)}"
        # The point of the fix: nothing was cut to a single character.
        assert all(len(code) > 1 for code in expected & codes)


@pytest.mark.python
@pytest.mark.storm
def test_hurdat_event_and_classification_stay_distinct():
    """A landfall event ``"L"`` and a low classification ``"LO"`` are different."""
    record = ("AL011990,          SYNTHTC,    2,\n"
              "19900801, 0000, L, LO,  25.0N,  80.0W,  30, 1005,"
              + "    0," * 12 + " -999,\n"
              "19900801, 0600,  , TS,  25.5N,  80.5W,  40, 1000,"
              + "    0," * 12 + " -999,\n")
    import tempfile
    from clawpack.geoclaw.met.track import iter_hurdat
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "codes.txt"
        path.write_text(record)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            track = next(iter_hurdat(path))

    assert track.event[0] == "L"
    assert track.classification[0] == "LO"
    assert track.classification[1] == "TS"


@pytest.mark.python
@pytest.mark.storm
def test_standardized_vocabularies():
    """Basin and status map onto enums across each format's own vocabulary."""
    from clawpack.geoclaw.met.track import (Basin, StormStatus,
                                            standardize_basin,
                                            standardize_status)

    # The same physical basin, spelled three ways by three archives.
    assert standardize_basin("AL", "atcf") is Basin.NORTH_ATLANTIC
    assert standardize_basin("AL", "hurdat") is Basin.NORTH_ATLANTIC
    assert standardize_basin("NA", "ibtracs") is Basin.NORTH_ATLANTIC
    assert standardize_basin("L", "tcvitals") is Basin.NORTH_ATLANTIC

    assert standardize_status("TS") is StormStatus.TROPICAL_STORM
    assert standardize_status("lo") is StormStatus.LOW
    # An unrecognized code degrades rather than stopping a read.
    assert standardize_basin("ZZ") is Basin.UNKNOWN
    assert standardize_status("ZZZ") is StormStatus.UNKNOWN
    assert standardize_basin(None) is Basin.UNKNOWN


@pytest.mark.python
@pytest.mark.storm
def test_track_exposes_standardized_accessors():
    """The string attributes are unchanged; the enums are reached separately."""
    pytest.importorskip("pandas")
    from clawpack.geoclaw.met.track import Basin, StormStatus

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        track = StormTrack.read_hurdat(_storm_input_path("hurdat"))

    # Source vocabulary preserved...
    assert track.basin == "Atlantic"
    assert track.basin_code == "AL"
    assert track.classification[0] == "LO"
    # ...and the canonical value available alongside it.
    assert track.basin_standard() is Basin.NORTH_ATLANTIC
    assert track.status()[0] is StormStatus.LOW


@pytest.mark.python
@pytest.mark.storm
def test_plot_swath_tolerates_missing_radius():
    """A track with no outer radius plots without a swath instead of raising.

    HURDAT2 supplies no storm radius at all, so ``plot(plot_swath=True)`` used
    to feed NaN straight into a matplotlib patch.
    """
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        track = StormTrack.read_hurdat(_storm_input_path("hurdat"))
    assert np.isnan(track.storm_radius).all()

    figure, ax = plt.subplots()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            track.plot(ax, plot_swath=True)
        assert any("swath is not drawn" in str(w.message) for w in caught)
        # An explicit radius still draws.
        track.plot(ax, plot_swath=True, radius=100e3)
    finally:
        plt.close(figure)
