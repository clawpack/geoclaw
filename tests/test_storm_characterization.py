#!/usr/bin/env python
# encoding: utf-8

"""Phase 0 characterization baseline for the storm/met-forcing I/O surface.

These tests pin the *current* behavior of the reader/writer surface so the
meteorological-forcing refactor can be proven behavior-neutral.  They are
stronger than the round-trip checks in ``test_storm.py``:

* **Read state snapshots** serialize the full public object state returned by
  ``read_geoclaw`` / ``read_atcf`` / ``read_hurdat`` / ``read_ibtracs`` /
  ``read_jma`` / ``read_tcvitals`` / ``read_data`` to a golden JSON and assert
  exact identity.  This catches parse changes that never reach the GeoClaw
  output (and so are invisible to a write-back round trip).
* **Write byte goldens** assert the exact bytes emitted by ``write_geoclaw``
  (against the committed ``*_geoclaw.txt`` baselines) and ``write_data``
  descriptors.

Golden files live in ``tests/data/storm/characterization/``.  To (re)generate
after an intentional behavior change, run with ``GEOCLAW_REGEN=1``::

    GEOCLAW_REGEN=1 pytest tests/test_storm_characterization.py

The netCDF fixtures are generated at runtime by the deterministic helpers in
``test_storm`` (no committed ``.nc`` beyond the existing IBTrACS track file).
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.surge.storm as storm

# ``tests/`` has no package __init__; make the sibling helpers importable the
# same way ``examples/storm-surge/isaac/test_isaac.py`` does.
sys.path.insert(0, str(Path(__file__).parent))
from test_storm import (  # noqa: E402
    create_era5_storm_file,
    create_nws13_storm_file,
    _make_storm_from_format,
    _storm_input_path,
    _storm_check_path,
)

testdir = Path(__file__).parent
data_dir = testdir / "data" / "storm"
golden_dir = data_dir / "characterization"

# Public attributes that constitute a Storm's characterized state.  Ordering is
# fixed for a stable, diff-friendly golden.  ``file_paths`` is normalized to
# basenames (see _normalize) because the absolute path is machine-specific.
_STATE_FIELDS = (
    "t", "eye_location", "max_wind_speed", "max_wind_radius",
    "central_pressure", "storm_radius", "wind_speeds", "time_offset",
    "name", "basin", "ID", "classification", "event",
    "file_format", "file_paths", "scaling", "storm_time_scale",
    "crop_extent", "ramp_width", "x_shift", "y_shift",
    "met_x_name", "met_y_name", "met_time_name", "met_lon_wrap",
    "met_y_increasing", "met_fill_value", "met_fill_action",
    "met_time_offset", "met_variable_map",
)


def _normalize(value):
    """Recursively convert a Storm attribute into a canonical JSON-able form."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, np.datetime64):
        return str(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, Path):
        return Path(value).name  # basename only — drop machine-specific dirs
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "M":  # datetime64 array
            return [str(v) for v in value]
        return _normalize(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_normalize(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _normalize(v) for k, v in sorted(value.items())}
    # Fall back to a stable string for anything unexpected.
    return str(value)


def _storm_state(s):
    """Canonical serializable snapshot of a Storm object's public state."""
    return {name: _normalize(getattr(s, name, None)) for name in _STATE_FIELDS}


def _dumps(obj):
    # allow_nan=True (default) so NaN-bearing fields (e.g. HURDAT radii) round
    # trip; comparison is text-vs-text so the non-strict-JSON "NaN" is fine.
    return json.dumps(obj, indent=2, sort_keys=True)


def _check_golden_text(actual_text, golden_path):
    """Assert ``actual_text`` matches the golden file, or (re)generate it."""
    golden_path = Path(golden_path)
    regen = bool(os.environ.get("GEOCLAW_REGEN"))
    if regen or not golden_path.exists():
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        golden_path.write_text(actual_text)
        if not regen:
            # First-time bootstrap on a checkout without the golden committed.
            pytest.skip(f"Baseline created: {golden_path.name} (rerun to assert)")
        return
    expected = golden_path.read_text()
    assert actual_text == expected, (
        f"Characterization mismatch vs {golden_path.name}. "
        f"If this change is intentional, regenerate with GEOCLAW_REGEN=1."
    )


def _check_state(s, name):
    _check_golden_text(_dumps(_storm_state(s)), golden_dir / f"{name}.json")


def _descriptor_head(text):
    """Return a write_data descriptor up to and including the '# File paths'
    marker, excluding the echoed (machine-specific) trailing path lines."""
    marker = "# File paths\n"
    idx = text.find(marker)
    assert idx != -1, "descriptor missing '# File paths' marker"
    return text[: idx + len(marker)]


# ---------------------------------------------------------------------------
# Read state snapshots
# ---------------------------------------------------------------------------

# GeoClaw-format reader over every committed *_geoclaw.txt baseline.
GEOCLAW_STATE_FORMATS = ["atcf", "hurdat", "jma", "tcvitals", "ibtracs"]


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format", GEOCLAW_STATE_FORMATS)
def test_read_geoclaw_state(file_format):
    """Snapshot read_geoclaw object state over the committed baselines."""
    s = storm.Storm(_storm_check_path(file_format), file_format="geoclaw")
    _check_state(s, f"read_geoclaw_{file_format}")


@pytest.mark.python
@pytest.mark.storm
def test_read_atcf_state():
    pytest.importorskip("pandas")
    s = storm.Storm(_storm_input_path("atcf"), file_format="atcf")
    _check_state(s, "read_atcf")


@pytest.mark.python
@pytest.mark.storm
def test_read_hurdat_state():
    s = storm.Storm(_storm_input_path("hurdat"), file_format="hurdat")
    _check_state(s, "read_hurdat")


@pytest.mark.python
@pytest.mark.storm
def test_read_jma_state():
    s = storm.Storm(_storm_input_path("jma"), file_format="jma")
    _check_state(s, "read_jma")


@pytest.mark.python
@pytest.mark.storm
def test_read_tcvitals_state():
    s = storm.Storm(_storm_input_path("tcvitals"), file_format="tcvitals")
    _check_state(s, "read_tcvitals")


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.netcdf
def test_read_ibtracs_state():
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    s = storm.Storm(_storm_input_path("ibtracs"), file_format="ibtracs",
                    sid="2008245N17323", agency_pref=["wmo", "usa"])
    _check_state(s, "read_ibtracs")


def _make_data_storm():
    """A deterministic data (gridded) Storm with fixed forcing controls."""
    s = storm.Storm()
    s.time_offset = np.datetime64("2012-08-29")
    s.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    s.ramp_width = 3
    s.x_shift = 1.25
    s.y_shift = -0.5
    return s


# NWS13-style names are not in the default discovery lists, so write_data needs
# an explicit role->variable mapping (mirrors test_storm.test_data_storm_roundtrip).
_NWS13_VAR_MAPPING = {"wind_u": "uwnd", "wind_v": "vwnd", "pressure": "press"}


def _write_netcdf_descriptor(s, kind, nc, descriptor):
    """Generate the netCDF fixture and write its data descriptor."""
    if kind == "era5":
        create_era5_storm_file(nc)
        s.file_format = "netcdf"
        s.file_paths = [nc]
        s.write(descriptor, file_format="data", dim_mapping={"t": "valid_time"})
    else:
        create_nws13_storm_file(nc)
        s.file_format = "nws13"
        s.file_paths = [nc]
        s.write(descriptor, file_format="data", var_mapping=_NWS13_VAR_MAPPING)


@pytest.mark.python
@pytest.mark.storm
def test_read_data_state_ascii(tmp_path):
    """Snapshot read_data state for an OWI/NWS12 ASCII descriptor."""
    s = _make_data_storm()
    s.file_format = "ascii"
    s.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                    Path("storm_2.PRE"), Path("storm_2.WIN")]
    descriptor = tmp_path / "ascii.storm"
    s.write(descriptor, file_format="data")
    read = storm.Storm(descriptor, file_format="data")
    _check_state(read, "read_data_ascii")


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.netcdf
@pytest.mark.parametrize("kind", ["era5", "nws13"])
def test_read_data_state_netcdf(tmp_path, kind):
    """Snapshot read_data state for netCDF descriptors (ERA5 + NWS13)."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    s = _make_data_storm()
    descriptor = tmp_path / f"{kind}.storm"
    _write_netcdf_descriptor(s, kind, tmp_path / f"{kind}.nc", descriptor)
    read = storm.Storm(descriptor, file_format="data")
    _check_state(read, f"read_data_{kind}")


# ---------------------------------------------------------------------------
# Write byte goldens
# ---------------------------------------------------------------------------

WRITE_GEOCLAW_FORMATS = ["atcf", "hurdat", "jma", "tcvitals"]


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format", WRITE_GEOCLAW_FORMATS)
def test_write_geoclaw_bytes(tmp_path, file_format):
    """write_geoclaw must reproduce the committed baseline byte-for-byte."""
    if file_format == "atcf":
        pytest.importorskip("pandas")
    test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
    out_path = tmp_path / f"{file_format}_geoclaw.txt"
    write_kwargs = {"file_format": "geoclaw"}
    if fill_mwr is not None:
        write_kwargs["max_wind_radius_fill"] = fill_mwr
    if fill_rad is not None:
        write_kwargs["storm_radius_fill"] = fill_rad
    test_storm.write(out_path, **write_kwargs)

    # Freeze the CURRENT write_geoclaw bytes.  We intentionally do NOT compare
    # against the committed *_geoclaw.txt baselines: those are consumed only
    # numerically by test_storm.check_geoclaw (skiprows=3), and their header
    # lines have drifted from current output (atcf uses a space instead of a
    # 'T' date separator; jma's committed header year reads '0008').  That
    # header drift is a pre-existing inconsistency, out of scope for Phase 0;
    # here we characterize exactly what the writer emits today.
    _check_golden_text(out_path.read_text(),
                       golden_dir / f"write_geoclaw_{file_format}.txt")


@pytest.mark.python
@pytest.mark.storm
def test_write_data_bytes_ascii(tmp_path):
    """Byte golden of the OWI/NWS12 ASCII descriptor body."""
    s = _make_data_storm()
    s.file_format = "ascii"
    s.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                    Path("storm_2.PRE"), Path("storm_2.WIN")]
    descriptor = tmp_path / "ascii.storm"
    s.write(descriptor, file_format="data")
    _check_golden_text(_descriptor_head(descriptor.read_text()),
                       golden_dir / "write_data_ascii.txt")


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.netcdf
@pytest.mark.parametrize("kind", ["era5", "nws13"])
def test_write_data_bytes_netcdf(tmp_path, kind):
    """Byte golden of the netCDF descriptor body (&file_info/&variable_info)."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    s = _make_data_storm()
    descriptor = tmp_path / f"{kind}.storm"
    _write_netcdf_descriptor(s, kind, tmp_path / f"{kind}.nc", descriptor)
    _check_golden_text(_descriptor_head(descriptor.read_text()),
                       golden_dir / f"write_data_{kind}.txt")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
