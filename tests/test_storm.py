#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing storm data."""

from pathlib import Path
import sys

import numpy as np
import pytest

import clawpack.geoclaw.surge.storm as storm

# Local test directory and bundled storm test data
testdir = Path(__file__).parent
data_dir = testdir / "data" / "storm"

# Current file-format tests
FILE_FORMAT_TESTS = ["atcf", 
                    "hurdat", 
                    "jma", 
                    "tcvitals", 
                    pytest.param("ibtracs", 
                                 marks=pytest.mark.xfail(
                                 reason=("IBTrACS GeoClaw output does not ,"
                                         "currently match baseline")))
]


def check_geoclaw(paths, check_header=False):
    """
    Check that two GeoClaw-formatted storm files are numerically equivalent.
    """
    paths = [Path(path) for path in paths]

    if check_header:
        with paths[0].open("r") as first_file, paths[1].open("r") as second_file:
            # Check for number of lines
            assert int(first_file.readline()) == int(second_file.readline())

            # Check for time offset
            assert first_file.readline() == second_file.readline()

    data = [np.loadtxt(path, skiprows=3) for path in paths]
    np.testing.assert_almost_equal(data[0], data[1])


def _storm_input_path(file_format):
    """Return the bundled input path for a given storm file format."""
    suffix = "nc" if file_format == "ibtracs" else "txt"
    return data_dir / f"{file_format}.{suffix}"


def _storm_check_path(file_format):
    """Return the bundled GeoClaw-format reference file for a storm format."""
    return data_dir / f"{file_format}_geoclaw.txt"


def _make_storm_from_format(file_format):
    """Read one bundled storm file and return the storm plus fill callbacks."""
    input_path = _storm_input_path(file_format)

    if file_format == "ibtracs":
        kwargs = {"sid": "2008245N17323", "agency_pref": ["wmo", "usa"]}

        atcf_path = data_dir / "atcf.txt"
        storm_atcf = storm.Storm(atcf_path, file_format="ATCF")

        def fill_mwr(t, this_storm):
            return storm.fill_rad_w_other_source(
                t, this_storm, storm_atcf, "max_wind_radius"
            )

        def fill_rad(t, this_storm):
            return storm.fill_rad_w_other_source(
                t, this_storm, storm_atcf, "storm_radius"
            )
    else:
        kwargs = {}
        fill_mwr = None
        fill_rad = None

    test_storm = storm.Storm(input_path, file_format=file_format, **kwargs)

    # Temporary normalization for formats that do not provide these radii.
    if file_format in ["hurdat", "jma"]:
        test_storm.max_wind_radius[:] = 0.0
        test_storm.storm_radius[:] = 0.0

    return test_storm, fill_mwr, fill_rad


@pytest.mark.python
@pytest.mark.parametrize("file_format", FILE_FORMAT_TESTS)
def test_storm_io(tmp_path, file_format):
    r"""Test reading one storm format and writing it back in GeoClaw format."""
    if file_format == "ibtracs":
        pytest.importorskip("xarray")
    if file_format == "atcf":
        pytest.importorskip("pandas")

    test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
    out_path = tmp_path / f"{file_format}_geoclaw.txt"
    check_path = _storm_check_path(file_format)
    write_kwargs = {"file_format": "geoclaw"}
    if fill_mwr is not None:
        write_kwargs["max_wind_radius_fill"] = fill_mwr
    if fill_rad is not None:
        write_kwargs["storm_radius_fill"] = fill_rad

    test_storm.write(out_path, **write_kwargs)

    check_geoclaw([out_path, check_path])


def save_storm_test_data(output_dir):
    """Utility helper to regenerate bundled GeoClaw-format storm baselines."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for file_format in FILE_FORMAT_TESTS:
        if file_format == "ibtracs":
            try:
                import xarray  # noqa: F401
            except ImportError:
                print("Skipping IBTrACS save, missing xarray.")
                continue
        if file_format == "atcf":
            try:
                import pandas  # noqa: F401
            except ImportError:
                print("Skipping ATCF save, missing pandas.")
                continue

        test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
        check_path = output_dir / f"{file_format}_geoclaw.txt"
        write_kwargs = {"file_format": "geoclaw"}
        if fill_mwr is not None:
            write_kwargs["max_wind_radius_fill"] = fill_mwr
        if fill_rad is not None:
            write_kwargs["storm_radius_fill"] = fill_rad
        test_storm.write(check_path, **write_kwargs)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "save":
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else data_dir
        save_storm_test_data(output_dir)
    else:
        pytest.main([__file__])
