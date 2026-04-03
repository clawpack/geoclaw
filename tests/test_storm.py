#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing storm data."""

from pathlib import Path
import sys

import numpy as np
import pytest

import clawpack.geoclaw.surge.storm as storm
import clawpack.geoclaw.util as util

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

DATA_FILE_FORMAT_MAP = {
    1: ["ascii", 1, "nws12", "owi"],
    2: ["netcdf", 2, "nws13"],
}


def create_netcdf_storm_file(path):
    """Create a small deterministic NetCDF storm dataset for testing."""
    xr = pytest.importorskip("xarray")
    size = (3, 5, 3)
    coords = [
        ("longitude", np.linspace(-1.0, 1.0, size[0])),
        ("latitude", np.linspace(-2.0, 2.0, size[1])),
        ("valid_time", np.linspace(0.0, 2.0, size[2])),
    ]

    values = np.arange(np.prod(size), dtype=float).reshape(size)
    wind_x = xr.DataArray(values, coords=coords)
    wind_y = xr.DataArray(values + 100.0, coords=coords)
    pressure = xr.DataArray(values + 1000.0, coords=coords)
    ds = xr.Dataset({"u": wind_x, "v": wind_y, "pressure": pressure})
    ds.to_netcdf(path)


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


@pytest.mark.python
@pytest.mark.parametrize("file_format", ["ascii", "netcdf"])
def test_data_storm_roundtrip(file_format, tmp_path):
    """Test round-trip read/write of data-driven storm metadata."""
    storm_path = tmp_path / "test.storm"
    data_storm = storm.Storm()
    data_storm.time_offset = np.datetime64("2012-08-29")
    data_storm.file_format = file_format
    data_storm.window_type = 1
    data_storm.ramp_width = 3

    if file_format == "ascii":
        data_storm.file_paths = [
            Path("storm_1.PRE"),
            Path("storm_1.WIN"),
            Path("storm_2.PRE"),
            Path("storm_2.WIN"),
        ]
        data_storm.write(storm_path, file_format="data")
        read_storm = storm.Storm(storm_path, file_format="data")
    else:
        data_storm.file_paths = [tmp_path / "storm.nc"]
        create_netcdf_storm_file(data_storm.file_paths[0])
        data_storm.write(
            storm_path,
            file_format="data",
            dim_mapping={"t": "valid_time"},
        )
        read_storm = storm.Storm(storm_path, file_format="data")

    assert data_storm.time_offset == read_storm.time_offset
    assert data_storm.file_format in DATA_FILE_FORMAT_MAP[read_storm.file_format]
    assert data_storm.window_type == read_storm.window_type
    assert data_storm.ramp_width == read_storm.ramp_width
    assert len(data_storm.file_paths) == len(read_storm.file_paths)
    for i, path in enumerate(data_storm.file_paths):
        assert read_storm.file_paths[i] == path


@pytest.mark.python
def test_netcdf_var_mapping(tmp_path):
    """Test NetCDF dimension and variable name discovery for data storms."""
    storm_data_file = tmp_path / "storm.nc"
    create_netcdf_storm_file(storm_data_file)

    dim_mapping = util.get_netcdf_names(
        storm_data_file,
        lookup_type="dim",
        verbose=True,
        user_mapping={"t": "valid_time"},
    )
    var_mapping = util.get_netcdf_names(
        storm_data_file,
        lookup_type="var",
        verbose=True,
    )

    assert dim_mapping == {"x": "longitude", "y": "latitude", "t": "valid_time"}
    assert var_mapping == {"wind_u": "u", "wind_v": "v", "pressure": "pressure"}


@pytest.mark.python
@pytest.mark.parametrize(
    "speeds_knots, expected_categories",
    [
        (np.array([20, 50, 70, 90, 100, 120, 140]), 
         np.array([-1,  0,  1,  2,   3,   4,   5])),
    ],
)
def test_storm_category_nhc(speeds_knots, expected_categories):
    """Test NHC categorization from known wind speeds."""
    s = storm.Storm()
    s.max_wind_speed = speeds_knots * 0.514444  # convert knots to m/s

    categories = s.category(categorization="NHC")

    print(f"{categories}")
    print(f"{expected_categories}")

    assert np.array_equal(categories, expected_categories)


@pytest.mark.python
@pytest.mark.parametrize(
    "speed_knots, expected_category, expected_name",
    [
        (20, -1, "Tropical Depression"),
        (50, 0, "Tropical Storm"),
        (80, 1, "Category 1 Hurricane"),
    ],
)
def test_storm_category_names(speed_knots, expected_category, expected_name):
    """Test NHC category-name output."""
    s = storm.Storm()
    s.max_wind_speed = np.array([speed_knots]) * 0.514444  # convert knots to m/s

    categories, names = s.category(categorization="NHC", cat_names=True)

    assert categories[0] == expected_category
    assert expected_name in names[0]


@pytest.mark.python
def test_storm_plot_smoke(tmp_path):
    """Smoke test for Storm.plot (no assertions, just ensure it runs)."""
    plt = pytest.importorskip("matplotlib.pyplot")

    s = storm.Storm()

    # Minimal valid track
    s.t = np.array([
        np.datetime64("2020-01-01"),
        np.datetime64("2020-01-02"),
    ])
    s.eye_location = np.array([
        [-90.0, 25.0],
        [-89.5, 25.5],
    ])

    fig, ax = plt.subplots()
    s.plot(ax)

    # Save to tmp_path just to ensure backend works
    fig.savefig(tmp_path / "storm_plot.png")
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "save":
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else data_dir
        save_storm_test_data(output_dir)
    else:
        raise SystemExit(pytest.main([__file__]))
