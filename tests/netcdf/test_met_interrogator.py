"""
Tests for MetInterrogator.

Covers: variable presence/absence, grid consistency, unit passthrough and
conversion, time offset calculation, ensemble dimension handling, and
descriptor output structure.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from clawpack.geoclaw.netcdf_utils import MetInterrogator, MetMetadata

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
    """All required variables present: interrogate_met() returns MetMetadata."""
    path = met_file_factory()
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()
    assert isinstance(meta, MetMetadata)


@pytest.mark.parametrize("missing_role", ["wind_u", "wind_v", "pressure"])
def test_missing_variable_raises(met_file_factory, missing_role):
    """
    interrogate_met() raises KeyError with a message identifying the missing
    variable when one required variable is absent from the file.
    """
    path = met_file_factory(omit_roles=[missing_role])
    expected_var = _VAR_MAP[missing_role]
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        with pytest.raises(KeyError, match=expected_var):
            intr.interrogate_met()


def test_empty_variable_map_raises(met_file_factory):
    """An empty variable_map raises ValueError before any file inspection."""
    path = met_file_factory()
    with MetInterrogator(path, variable_map={}) as intr:
        with pytest.raises(ValueError, match="variable_map"):
            intr.interrogate_met()


# ============================================================
# Grid / dimension consistency
# ============================================================

def test_mismatched_variable_dims_raises(met_file_factory):
    """
    interrogate_met() raises ValueError when variables do not share
    the same set of dimensions.
    """
    path = met_file_factory(mismatch_dims=True)
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        with pytest.raises(ValueError, match="[Gg]rid|[Dd]im"):
            intr.interrogate_met()


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
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        with pytest.raises(ValueError, match="[Tt]ime"):
            intr.interrogate_met()


# ============================================================
# Ensemble dimension handling
# ============================================================

def test_singleton_ensemble_dim_passes(met_file_factory):
    """
    A singleton extra dimension (e.g. a single ensemble member) does not
    cause interrogate_met() to raise.
    """
    path = met_file_factory(extra_dim=("member", 1))
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()
    assert isinstance(meta, MetMetadata)


def test_non_singleton_ensemble_dim_raises(met_file_factory):
    """
    A non-singleton extra dimension raises ValueError — GeoClaw does not
    support ensemble/member dimensions.
    """
    path = met_file_factory(extra_dim=("member", 3))
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        with pytest.raises(ValueError, match="[Ee]nsemble|[Mm]ember|singleton"):
            intr.interrogate_met()


# ============================================================
# Unit verification and passthrough
# ============================================================

@pytest.mark.parametrize("unit_cfg", PRESSURE_UNIT_VARIANTS)
def test_pressure_unit_variants(met_file_factory, unit_cfg):
    """
    Pressure unit recorded in MetVariableInfo.source_units matches the unit
    string in the file.  Conversion-needed cases emit a warning; Pa passes
    silently.
    """
    pressure_units = unit_cfg["pressure_units"]
    expected_source = unit_cfg["expected_source"]

    path = met_file_factory(pressure_units=pressure_units)
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()

    pressure_info = next(v for v in meta.variables if v.geoclaw_role == "pressure")
    assert pressure_info.source_units == expected_source, (
        f"For pressure_units={pressure_units!r}: expected source_units="
        f"{expected_source!r}, got {pressure_info.source_units!r}"
    )


@pytest.mark.parametrize("unit_cfg", WIND_UNIT_VARIANTS)
def test_wind_unit_variants(met_file_factory, unit_cfg):
    """
    Wind unit recorded in MetVariableInfo.source_units matches the unit
    string in the file.  knots triggers a warning; m/s passes silently.
    """
    wind_units = unit_cfg["wind_units"]
    expected_source = unit_cfg["expected_source"]

    path = met_file_factory(wind_units=wind_units)
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()

    wind_u_info = next(v for v in meta.variables if v.geoclaw_role == "wind_u")
    assert wind_u_info.source_units == expected_source


def test_contract_units_no_warning(met_file_factory):
    """
    Files already in contract units (Pa, m/s) do not emit any unit warnings.
    """
    path = met_file_factory(pressure_units="Pa", wind_units="m/s")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
            intr.interrogate_met()


def test_unrecognised_pressure_unit_raises(met_file_factory):
    """An unrecognised pressure unit string raises ValueError."""
    path = met_file_factory(pressure_units="atm_weird_unknown")
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        with pytest.raises(ValueError, match="[Uu]nrecogni"):
            intr.interrogate_met()


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
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()

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
    with MetInterrogator(path, variable_map=_VAR_MAP,
                          time_reference=time_start) as intr:
        meta = intr.interrogate_met()
    assert meta.time_offset == pytest.approx(0.0, abs=1.0), (
        f"Expected time_offset≈0, got {meta.time_offset}"
    )


# ============================================================
# MetMetadata structure
# ============================================================

def test_metadata_variables_all_roles_present(met_file_factory):
    """
    MetMetadata.variables contains one entry per variable_map role.
    """
    path = met_file_factory()
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()

    roles = {v.geoclaw_role for v in meta.variables}
    assert roles == set(_VAR_MAP.keys()), (
        f"Expected roles {set(_VAR_MAP.keys())}, got {roles}"
    )


def test_metadata_variable_names_match_map(met_file_factory):
    """MetVariableInfo.var_name matches the value in variable_map."""
    path = met_file_factory()
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()

    for v in meta.variables:
        assert v.var_name == _VAR_MAP[v.geoclaw_role], (
            f"Role {v.geoclaw_role!r}: expected var_name="
            f"{_VAR_MAP[v.geoclaw_role]!r}, got {v.var_name!r}"
        )


def test_fill_action_default_is_warn(met_file_factory):
    """Default fill_action for MetInterrogator is 'warn'."""
    path = met_file_factory()
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()
    assert meta.fill_action == "warn"


def test_fill_action_abort_propagated(met_file_factory):
    """fill_action='abort' is forwarded into MetMetadata."""
    path = met_file_factory()
    with MetInterrogator(path, variable_map=_VAR_MAP,
                          fill_action="abort") as intr:
        meta = intr.interrogate_met()
    assert meta.fill_action == "abort"


def test_invalid_fill_action_raises(met_file_factory):
    """
    Passing an invalid fill_action raises ValueError at construction.

    A real file path is required because MetInterrogator opens the file in
    super().__init__() before the fill_action validation check fires.
    """
    path = met_file_factory()
    with pytest.raises(ValueError, match="fill_action"):
        MetInterrogator(path, variable_map=_VAR_MAP, fill_action="ignore")


# ============================================================
# All coordinate and dim-order variants
# ============================================================

@pytest.mark.parametrize("coord_kwargs", COORD_VARIANTS)
def test_coord_variants(met_file_factory, coord_kwargs):
    """interrogate_met() completes for every coordinate variant."""
    path = met_file_factory(**coord_kwargs)
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()
    expected_conv = 360 if coord_kwargs["lon_max"] > 180 else 180
    assert meta.lon_convention == expected_conv


@pytest.mark.parametrize("dim_kwargs", MET_DIM_ORDER_VARIANTS)
def test_dim_order_variants(met_file_factory, dim_kwargs):
    """interrogate_met() reports correct dim_order for all axis arrangements."""
    path = met_file_factory(**dim_kwargs)
    with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
        meta = intr.interrogate_met()
    expected = [{"lat": "lat", "lon": "lon", "time": "time"}.get(d, d)
                for d in dim_kwargs["dim_order"]]
    assert meta.dim_order == expected
