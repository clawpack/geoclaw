"""
Tests for DescriptorWriter.

Covers: topo descriptor round-trip, met descriptor structure, optional
fields omitted when not specified, and crop_bounds presence/absence.
"""
import io
import re

import pytest

from clawpack.geoclaw.netcdf_utils import (
    DescriptorWriter,
    TopoInterrogator,
    MetInterrogator,
    TopoMetadata,
    MetMetadata,
    MetVariableInfo,
)

from ._helpers import make_topo_dataset, make_met_dataset

pytestmark = [pytest.mark.python, pytest.mark.netcdf]

_VAR_MAP = {"wind_u": "u10", "wind_v": "v10", "pressure": "msl"}


# ============================================================
# Helpers
# ============================================================

def _parse_topo_descriptor(text: str) -> dict:
    """
    Parse a topo descriptor block (key = value lines) into a dict.
    Stops at the first blank line.
    """
    result = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            break
        if "=" in line:
            key, _, value = line.partition("=")
            result[key.strip()] = value.strip()
    return result


def _parse_met_descriptor(text: str) -> dict:
    """
    Parse a met descriptor written by write_met_descriptor into:
      {
        'file_info': {key: value, ...},
        'variable_info': [ {key: value, ...}, ... ]
      }
    """
    result: dict = {"file_info": {}, "variable_info": []}
    in_file_info = False
    in_var_info = False
    current_var: dict = {}

    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "&file_info":
            in_file_info = True
            in_var_info = False
            continue
        if stripped.startswith("&variable_info"):
            if in_var_info and current_var:
                result["variable_info"].append(current_var)
            in_file_info = False
            in_var_info = True
            current_var = {}
            # variable_info content may be on the same line
            rest = stripped[len("&variable_info"):].strip()
            for kv in re.findall(r'(\w+)=(\S+)', rest):
                current_var[kv[0]] = kv[1]
            continue
        if stripped == "/":
            if in_file_info:
                in_file_info = False
            elif in_var_info:
                if current_var:
                    result["variable_info"].append(current_var)
                current_var = {}
                in_var_info = False
            continue
        if in_file_info and "=" in stripped:
            key, _, value = stripped.partition("=")
            result["file_info"][key.strip()] = value.strip()
        if in_var_info and "=" in stripped:
            for kv in re.findall(r'(\w+)=(\S+)', stripped):
                current_var[kv[0]] = kv[1]

    if in_var_info and current_var:
        result["variable_info"].append(current_var)

    return result


def _get_topo_meta(topo_file_factory, **kwargs) -> TopoMetadata:
    """Build a TopoMetadata via interrogation for use in writer tests."""
    import warnings
    path = topo_file_factory(**kwargs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with TopoInterrogator(path, var_name="z") as intr:
            return intr.interrogate_topo()


def _get_met_meta(met_file_factory, **kwargs) -> MetMetadata:
    """Build a MetMetadata via interrogation for use in writer tests."""
    import warnings
    path = met_file_factory(**kwargs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with MetInterrogator(path, variable_map=_VAR_MAP) as intr:
            return intr.interrogate_met()


# ============================================================
# Topo descriptor
# ============================================================

def test_topo_descriptor_required_keys_present(topo_file_factory):
    """
    write_topo_descriptor emits all mandatory keys: var_name, lon_name,
    lat_name, lon_offset, lat_order, dim_order, fill_action.
    """
    meta = _get_topo_meta(topo_file_factory)
    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    parsed = _parse_topo_descriptor(buf.getvalue())

    for key in ("var_name", "lon_name", "lat_name", "lon_offset",
                "lat_order", "dim_order", "fill_action"):
        assert key in parsed, f"Missing required key '{key}' in topo descriptor"


def test_topo_descriptor_values_match_metadata(topo_file_factory):
    """
    Values in the topo descriptor text match the TopoMetadata they came from.
    """
    meta = _get_topo_meta(topo_file_factory,
                          lon_min=-10.0, lon_max=10.0,
                          lat_min=-10.0, lat_max=10.0,
                          lat_direction="N_to_S")
    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    parsed = _parse_topo_descriptor(buf.getvalue())

    assert parsed["var_name"] == meta.var_name
    assert parsed["lon_name"] == meta.lon_name
    assert parsed["lat_name"] == meta.lat_name
    assert float(parsed["lon_offset"]) == pytest.approx(meta.lon_offset)
    assert parsed["lat_order"] == meta.lat_order
    assert parsed["dim_order"] == ",".join(meta.dim_order)
    assert parsed["fill_action"] == meta.fill_action


def test_topo_descriptor_fill_value_omitted_when_none(topo_file_factory):
    """fill_value line is absent from descriptor when TopoMetadata.fill_value is None."""
    meta = _get_topo_meta(topo_file_factory)
    # Force fill_value to None in the metadata
    import dataclasses
    meta = dataclasses.replace(meta, fill_value=None)

    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    parsed = _parse_topo_descriptor(buf.getvalue())

    assert "fill_value" not in parsed


def test_topo_descriptor_fill_value_present_when_set(topo_file_factory):
    """fill_value line IS written when TopoMetadata.fill_value is not None."""
    import dataclasses
    meta = _get_topo_meta(topo_file_factory)
    meta = dataclasses.replace(meta, fill_value=-9999.0)

    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    parsed = _parse_topo_descriptor(buf.getvalue())

    assert "fill_value" in parsed
    assert float(parsed["fill_value"]) == pytest.approx(-9999.0)


def test_topo_descriptor_crop_bounds_omitted_when_none(topo_file_factory):
    """crop_bounds line is absent when no crop_bounds were specified."""
    meta = _get_topo_meta(topo_file_factory)
    assert meta.crop_bounds is None   # confirm fixture has no crop

    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    parsed = _parse_topo_descriptor(buf.getvalue())

    assert "crop_bounds" not in parsed


def test_topo_descriptor_crop_bounds_present_when_set(topo_file_factory):
    """crop_bounds line is written with all four values when specified."""
    import dataclasses
    meta = _get_topo_meta(topo_file_factory,
                          lon_min=-10.0, lon_max=10.0,
                          lat_min=-10.0, lat_max=10.0)
    meta = dataclasses.replace(meta, crop_bounds=(-5.0, 5.0, -5.0, 5.0))

    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    parsed = _parse_topo_descriptor(buf.getvalue())

    assert "crop_bounds" in parsed
    parts = parsed["crop_bounds"].split()
    assert len(parts) == 4
    lon0, lon1, lat0, lat1 = map(float, parts)
    assert lon0 == pytest.approx(-5.0)
    assert lon1 == pytest.approx(5.0)
    assert lat0 == pytest.approx(-5.0)
    assert lat1 == pytest.approx(5.0)


def test_descriptor_writes_lon_offset_not_convention(topo_file_factory):
    """
    Topo descriptor must contain 'lon_offset' as a parseable float and must
    NOT contain 'lon_convention'.  This verifies the format-version-2 change.
    """
    import dataclasses
    meta = _get_topo_meta(topo_file_factory)
    # Exercise a non-zero offset to confirm the value round-trips.
    meta = dataclasses.replace(meta, lon_offset=-360.0)

    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    text = buf.getvalue()
    parsed = _parse_topo_descriptor(text)

    assert "lon_offset" in parsed, "Descriptor must contain 'lon_offset'"
    assert "lon_convention" not in parsed, (
        "Descriptor must NOT contain 'lon_convention' (superseded by lon_offset)"
    )
    assert float(parsed["lon_offset"]) == pytest.approx(-360.0)


def test_topo_descriptor_ends_with_blank_line(topo_file_factory):
    """
    Descriptor block ends with a blank line so the Fortran parser can detect
    the end of the block.
    """
    meta = _get_topo_meta(topo_file_factory)
    buf = io.StringIO()
    DescriptorWriter.write_topo_descriptor(buf, meta)
    text = buf.getvalue()
    assert text.endswith("\n\n"), (
        "Topo descriptor must end with a blank line (two consecutive newlines) "
        "for the Fortran read loop to terminate correctly"
    )


# ============================================================
# Met descriptor
# ============================================================

def test_met_descriptor_file_info_block_present(met_file_factory):
    """write_met_descriptor writes a &file_info block."""
    meta = _get_met_meta(met_file_factory)
    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())

    assert parsed["file_info"], "file_info block is missing or empty"


def test_met_descriptor_file_info_required_keys(met_file_factory):
    """
    &file_info block contains all mandatory keys: lon_name, lat_name,
    time_name, dim_order, lon_convention, lat_order, fill_action, time_offset.
    """
    meta = _get_met_meta(met_file_factory)
    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())
    fi = parsed["file_info"]

    for key in ("lon_name", "lat_name", "time_name",
                "dim_order", "lon_convention", "lat_order",
                "fill_action", "time_offset"):
        assert key in fi, f"Missing key '{key}' in &file_info block"


def test_met_descriptor_variable_info_blocks_for_all_roles(met_file_factory):
    """
    One &variable_info block is written for each entry in variable_map.
    """
    meta = _get_met_meta(met_file_factory)
    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())

    assert len(parsed["variable_info"]) == len(_VAR_MAP), (
        f"Expected {len(_VAR_MAP)} variable_info blocks, "
        f"got {len(parsed['variable_info'])}"
    )


def test_met_descriptor_variable_info_roles_present(met_file_factory):
    """
    Each &variable_info block contains geoclaw_role and var_name.
    """
    meta = _get_met_meta(met_file_factory)
    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())

    for block in parsed["variable_info"]:
        assert "var_name" in block, (
            f"variable_info block missing 'var_name': {block}"
        )
        assert "geoclaw_role" in block, (
            f"variable_info block missing 'geoclaw_role': {block}"
        )


def test_met_descriptor_crop_bounds_omitted_when_none(met_file_factory):
    """crop_bounds is absent from &file_info when not specified."""
    meta = _get_met_meta(met_file_factory)
    assert meta.crop_bounds is None

    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())

    assert "crop_bounds" not in parsed["file_info"]


def test_met_descriptor_crop_bounds_present_when_set(met_file_factory):
    """crop_bounds appears in &file_info when crop_bounds are specified."""
    import dataclasses
    meta = _get_met_meta(met_file_factory,
                         lon_min=-10.0, lon_max=10.0,
                         lat_min=-10.0, lat_max=10.0)
    meta = dataclasses.replace(meta, crop_bounds=(-5.0, 5.0, -5.0, 5.0))

    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())

    assert "crop_bounds" in parsed["file_info"]


def test_met_descriptor_fill_value_omitted_when_none(met_file_factory):
    """fill_value line is absent from &file_info when MetMetadata.fill_value is None."""
    import dataclasses
    meta = _get_met_meta(met_file_factory)
    meta = dataclasses.replace(meta, fill_value=None)

    buf = io.StringIO()
    DescriptorWriter.write_met_descriptor(buf, meta)
    parsed = _parse_met_descriptor(buf.getvalue())

    assert "fill_value" not in parsed["file_info"]
