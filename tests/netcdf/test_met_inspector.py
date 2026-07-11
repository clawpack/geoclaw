"""
Tests for MetInspector.

Covers: variable presence/absence, grid consistency, unit passthrough and
conversion, time offset calculation, ensemble dimension handling, and
descriptor output structure.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from clawpack.geoclaw.netcdf_utils import MetInspector, MetMetadata

from ._helpers import (
    make_met_dataset,
    COORD_VARIANTS,
    MET_DIM_ORDER_VARIANTS,
    PRESSURE_UNIT_VARIANTS,
    WIND_UNIT_VARIANTS,
)

pytestmark = [pytest.mark.python, pytest.mark.netcdf]

# Default variable map used throughout
_VAR_MAP = {"wind_u": "u10", "wind_v": "v10", "pressure": "msl"}


# ============================================================
# Variable presence / absence
# ============================================================

def test_all_variables_present_passes(met_file_factory):
    """All required variables present: inspect_met() returns MetMetadata."""
    path = met_file_factory()
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()
    assert isinstance(meta, MetMetadata)


@pytest.mark.parametrize("missing_role", ["wind_u", "wind_v", "pressure"])
def test_missing_variable_raises(met_file_factory, missing_role):
    """
    inspect_met() raises ValueError identifying the missing variable when an
    explicit variable_map points to a name absent from the file.
    """
    path = met_file_factory(omit_roles=[missing_role])
    expected_var = _VAR_MAP[missing_role]
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match=expected_var):
            insp.inspect_met()


def test_empty_variable_map_autodiscovers(met_file_factory):
    """An empty (or omitted) variable_map auto-discovers all met roles.

    The default met file uses common variable names (u10/v10/msl) that are
    resolved from the fallback lists, so no explicit map is needed.
    """
    path = met_file_factory()
    with MetInspector(path, variable_map={}) as insp:
        meta = insp.inspect_met()
    roles = {v.geoclaw_role: v.var_name for v in meta.variables}
    assert roles == {"wind_u": "u10", "wind_v": "v10", "pressure": "msl"}


# ============================================================
# Grid / dimension consistency
# ============================================================

def test_mismatched_variable_dims_raises(met_file_factory):
    """
    inspect_met() raises ValueError when variables do not share
    the same set of dimensions.
    """
    path = met_file_factory(mismatch_dims=True)
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="[Gg]rid|[Dd]im"):
            insp.inspect_met()


def test_missing_time_coord_raises(met_file_factory):
    """
    A file with no time coordinate raises ValueError (met forcing requires
    a time axis).
    """
    import xarray as xr
    ds = make_met_dataset()
    # Remove the time dimension by selecting a single time step
    ds_no_time = ds.isel(time=0).drop_vars("time")
    path = met_file_factory(ds=ds_no_time)
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="[Tt]ime"):
            insp.inspect_met()


# ============================================================
# Ensemble dimension handling
# ============================================================

def test_singleton_ensemble_dim_passes(met_file_factory):
    """
    A singleton extra dimension (e.g. a single ensemble member) does not
    cause inspect_met() to raise.
    """
    path = met_file_factory(extra_dim=("member", 1))
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()
    assert isinstance(meta, MetMetadata)


def test_non_singleton_ensemble_dim_raises(met_file_factory):
    """
    A non-singleton extra dimension raises ValueError — GeoClaw does not
    support ensemble/member dimensions.
    """
    path = met_file_factory(extra_dim=("member", 3))
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="[Ee]nsemble|[Mm]ember|singleton"):
            insp.inspect_met()


# ============================================================
# Unit verification and passthrough
# ============================================================

@pytest.mark.parametrize("unit_cfg", PRESSURE_UNIT_VARIANTS)
def test_pressure_unit_variants(met_file_factory, unit_cfg):
    """
    Pressure in the contract unit (Pa) passes and is recorded in
    source_units; any non-contract unit (hPa, mbar) is rejected rather than
    silently misread, since the read path does not convert.
    """
    pressure_units = unit_cfg["pressure_units"]
    expected_source = unit_cfg["expected_source"]

    path = met_file_factory(pressure_units=pressure_units)
    if pressure_units == "Pa":
        with MetInspector(path, variable_map=_VAR_MAP) as insp:
            meta = insp.inspect_met()
        pressure_info = next(v for v in meta.variables
                             if v.geoclaw_role == "pressure")
        assert pressure_info.source_units == expected_source
    else:
        with MetInspector(path, variable_map=_VAR_MAP) as insp:
            with pytest.raises(ValueError, match=pressure_units):
                insp.inspect_met()


@pytest.mark.parametrize("unit_cfg", WIND_UNIT_VARIANTS)
def test_wind_unit_variants(met_file_factory, unit_cfg):
    """
    Wind in the contract unit (m/s) passes and is recorded in source_units;
    a non-contract unit (knots) is rejected rather than silently misread,
    since the read path does not convert.
    """
    wind_units = unit_cfg["wind_units"]
    expected_source = unit_cfg["expected_source"]

    path = met_file_factory(wind_units=wind_units)
    if wind_units == "m/s":
        with MetInspector(path, variable_map=_VAR_MAP) as insp:
            meta = insp.inspect_met()
        wind_u_info = next(v for v in meta.variables
                           if v.geoclaw_role == "wind_u")
        assert wind_u_info.source_units == expected_source
    else:
        with MetInspector(path, variable_map=_VAR_MAP) as insp:
            with pytest.raises(ValueError, match="m/s|" + wind_units):
                insp.inspect_met()


def test_contract_units_no_warning(met_file_factory):
    """
    Files already in contract units (Pa, m/s) do not emit any unit warnings.
    """
    path = met_file_factory(pressure_units="Pa", wind_units="m/s")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with MetInspector(path, variable_map=_VAR_MAP) as insp:
            insp.inspect_met()


def test_unrecognised_pressure_unit_raises(met_file_factory):
    """An unrecognised pressure unit string raises ValueError."""
    path = met_file_factory(pressure_units="atm_weird_unknown")
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="[Uu]nrecogni"):
            insp.inspect_met()


def test_missing_units_raises(met_file_factory):
    """A forcing variable with no units attribute raises -- units are required
    and never silently assumed."""
    path = met_file_factory(wind_units="", pressure_units="")
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="no 'units'"):
            insp.inspect_met()


def test_missing_units_assume_units_opt_in(met_file_factory):
    """The explicit assume_units escape hatch assumes each variable's contract
    unit for a unitless file, without raising or warning."""
    path = met_file_factory(wind_units="", pressure_units="")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with MetInspector(path, variable_map=_VAR_MAP,
                          assume_units=True) as insp:
            meta = insp.inspect_met()
    roles = {v.geoclaw_role: v.source_units for v in meta.variables}
    assert roles["pressure"] == "Pa"
    assert roles["wind_u"] == "m/s"


# ============================================================
# Time offset calculation
# ============================================================

def test_time_offset_unix_epoch_default(met_file_factory):
    """
    With no time_reference, time_offset is seconds between the Unix epoch
    (1970-01-01) and the first time value in the file.
    """
    time_start = "2020-01-01"
    path = met_file_factory(time_start=time_start)
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()

    expected_offset = pd.Timestamp(time_start).timestamp()
    assert meta.time_offset == pytest.approx(expected_offset, rel=1e-6), (
        f"Expected time_offset≈{expected_offset}, got {meta.time_offset}"
    )


def test_time_offset_custom_reference(met_file_factory):
    """
    time_offset is zero when time_reference matches the first time in the file.
    """
    time_start = "2020-06-15T12:00:00"
    path = met_file_factory(time_start=time_start)
    with MetInspector(path, variable_map=_VAR_MAP,
                          time_reference=time_start) as insp:
        meta = insp.inspect_met()
    assert meta.time_offset == pytest.approx(0.0, abs=1.0), (
        f"Expected time_offset≈0, got {meta.time_offset}"
    )


def test_numeric_time_axis_raises(tmp_path):
    """A bare numeric met time axis (a duration with no 'since' reference) is
    rejected rather than silently interpreted as nanoseconds by pd.Timestamp."""
    ds = make_met_dataset()
    ds = ds.assign_coords(time=np.array([0.0, 6.0, 12.0]))
    ds["time"].attrs["units"] = "hours"
    ds["time"].encoding = {}
    path = tmp_path / "met_numeric_time.nc"
    ds.to_netcdf(path)
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="absolute time reference"):
            insp.inspect_met()


def test_non_seconds_time_axis_raises(tmp_path):
    """The met time axis must be 'seconds since <date>'.  An 'hours since'
    (or other non-second) datetime axis is rejected, because the Fortran met
    reader treats raw time values as integer seconds without unit scaling."""
    ds = make_met_dataset()
    ds["time"].encoding = {}
    path = tmp_path / "met_hours.nc"
    ds.to_netcdf(path, encoding={"time": {"units": "hours since 2020-01-01"}})
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        with pytest.raises(ValueError, match="requires the time axis"):
            insp.inspect_met()


# ============================================================
# MetMetadata structure
# ============================================================

def test_metadata_variables_all_roles_present(met_file_factory):
    """
    MetMetadata.variables contains one entry per variable_map role.
    """
    path = met_file_factory()
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()

    roles = {v.geoclaw_role for v in meta.variables}
    assert roles == set(_VAR_MAP.keys()), (
        f"Expected roles {set(_VAR_MAP.keys())}, got {roles}"
    )


def test_metadata_variable_names_match_map(met_file_factory):
    """MetVariableInfo.var_name matches the value in variable_map."""
    path = met_file_factory()
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()

    for v in meta.variables:
        assert v.var_name == _VAR_MAP[v.geoclaw_role], (
            f"Role {v.geoclaw_role!r}: expected var_name="
            f"{_VAR_MAP[v.geoclaw_role]!r}, got {v.var_name!r}"
        )


def test_fill_action_default_is_warn(met_file_factory):
    """Default fill_action for MetInspector is 'warn'."""
    path = met_file_factory()
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()
    assert meta.fill_action == "warn"


def test_fill_action_abort_propagated(met_file_factory):
    """fill_action='abort' is forwarded into MetMetadata."""
    path = met_file_factory()
    with MetInspector(path, variable_map=_VAR_MAP,
                          fill_action="abort") as insp:
        meta = insp.inspect_met()
    assert meta.fill_action == "abort"


def test_invalid_fill_action_raises(met_file_factory):
    """
    Passing an invalid fill_action raises ValueError at construction.

    A real file path is required because MetInspector opens the file in
    super().__init__() before the fill_action validation check fires.
    """
    path = met_file_factory()
    with pytest.raises(ValueError, match="fill_action"):
        MetInspector(path, variable_map=_VAR_MAP, fill_action="ignore")


# ============================================================
# All coordinate and dim-order variants
# ============================================================

@pytest.mark.parametrize("coord_kwargs", COORD_VARIANTS)
def test_coord_variants(met_file_factory, coord_kwargs):
    """inspect_met() completes for every coordinate variant."""
    path = met_file_factory(**coord_kwargs)
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()
    expected_conv = 360 if coord_kwargs["lon_max"] > 180 else 180
    assert meta.lon_wrap == expected_conv


@pytest.mark.parametrize("dim_kwargs", MET_DIM_ORDER_VARIANTS)
def test_dim_order_variants(met_file_factory, dim_kwargs):
    """inspect_met() reports correct dim_order for all axis arrangements."""
    path = met_file_factory(**dim_kwargs)
    with MetInspector(path, variable_map=_VAR_MAP) as insp:
        meta = insp.inspect_met()
    expected = [{"lat": "y", "lon": "x", "time": "time"}.get(d, d)
                for d in dim_kwargs["dim_order"]]
    assert meta.dim_order == expected
