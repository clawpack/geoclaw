#!/usr/bin/env python

import sys
from pathlib import Path
from urllib.error import URLError

import numpy as np
import pytest

import clawpack.clawutil.data
import clawpack.geoclaw.topotools as topotools

# Local test directory and bundled test data
testdir = Path(__file__).parent
data_dir = testdir / "data"
etopo1_extent = [-125.0, -124.0, 48.0, 48.5]

# Test topography functions
def topo_bowl(x, y):
    """Sample bowl topography."""
    return 1000.0 * (x**2 + y**2 - 1.0)


def topo_bowl_hill(x, y):
    """Sample bowl topography with a Gaussian hill."""
    z = 1000.0 * (x**2 + y**2 - 1.0)
    z = z + 1000.0 * np.exp(-100.0 * ((x - 0.7) ** 2 + (y - 0.8) ** 2))
    return z


@pytest.mark.python
@pytest.mark.parametrize("topo_type", [1, 2, 3])
def test_read_write_topo_bowl(tmp_path, topo_type):
    """
    Test writing and reading topo files with a small number of points.
    Note that ordering should go from the NW corner.
    """

    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = np.linspace(-1.0, 3.0, 5)
    topo.y = np.linspace(0.0, 3.0, 4)

    assert np.allclose(topo.x, np.array([-1.0, 0.0, 1.0, 2.0, 3.0])), \
        "Topography x values are incorrect."
    assert np.allclose(
        topo.X,
        np.array(
            [
                [-1.0, 0.0, 1.0, 2.0, 3.0],
                [-1.0, 0.0, 1.0, 2.0, 3.0],
                [-1.0, 0.0, 1.0, 2.0, 3.0],
                [-1.0, 0.0, 1.0, 2.0, 3.0],
            ]
        ),
    ), "Topography X values are incorrect."
    assert np.allclose(topo.y, np.array([0.0, 1.0, 2.0, 3.0])), \
        "Topography y values are incorrect."
    assert np.allclose(
        topo.Y,
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0, 2.0, 2.0],
                [3.0, 3.0, 3.0, 3.0, 3.0],
            ]
        ),
    ), "Topography Y values are incorrect."
    assert np.allclose(
        topo.Z,
        np.array(
            [
                [0.0, -1000.0, 0.0, 3000.0, 8000.0],
                [1000.0, 0.0, 1000.0, 4000.0, 9000.0],
                [4000.0, 3000.0, 4000.0, 7000.0, 12000.0],
                [9000.0, 8000.0, 9000.0, 12000.0, 17000.0],
            ]
        ),
    ), "Topography Z values are incorrect."

    path = tmp_path / f"bowl.tt{topo_type}"
    topo.write(path, topo_type=topo_type, Z_format="%22.15e")

    topo_in = topotools.Topography(path)
    assert np.allclose(topo.Z, topo_in.Z), \
        "Difference in written and read topography found."


@pytest.mark.python
def test_crop_topo_bowl():
    """
    Test cropping a topo file.
    """

    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = np.linspace(-1.0, 3.0, 5)
    topo.y = np.linspace(0.0, 3.0, 4)

    # topo.Z should be created automatically when referenced below:
    assert np.allclose(
        topo.Z,
        np.array(
            [
                [0.0, -1000.0, 0.0, 3000.0, 8000.0],
                [1000.0, 0.0, 1000.0, 4000.0, 9000.0],
                [4000.0, 3000.0, 4000.0, 7000.0, 12000.0],
                [9000.0, 8000.0, 9000.0, 12000.0, 17000.0],
            ]
        ),
    ), "Basic topography does not match test data."

    cropped_topo = topo.crop([0, 1, 0, 2])
    assert np.allclose(cropped_topo.x, np.array([0.0, 1.0])), \
        "Cropped topography y values do not match"
    assert np.allclose(cropped_topo.y, np.array([0.0, 1.0, 2.0])), \
        "Cropped topography y values do not match."
    assert np.allclose(
        cropped_topo.Z,
        np.array(
            [
                [-1000.0, 0.0],
                [0.0, 1000.0],
                [3000.0, 4000.0],
            ]
        ),
    ), "Cropped topography Z values do not match."


@pytest.mark.python
@pytest.mark.parametrize("topo_type", [1, 2, 3])
def test_read_write_topo_bowl_hill(tmp_path, topo_type):
    """Test writing and reading topo files for the bowl-hill example."""

    topo = topotools.Topography(topo_func=topo_bowl_hill)
    topo.x = np.linspace(-1.5, 2.5, 101)
    topo.y = np.linspace(-1.0, 2.0, 76)

    file_path = tmp_path / f"bowl_hill.tt{topo_type}"
    topo.write(file_path, topo_type=topo_type, Z_format="%22.15e")
    topo_in = topotools.Topography(path=file_path, topo_type=topo_type)

    assert np.allclose(topo.Z, topo_in.Z), (
        f"Written file of topo_type={topo_type} does not equal read in file."
    )


@pytest.mark.python
@pytest.mark.netcdf
def test_netcdf(tmp_path):
    r"""Test NetCDF topography I/O using a checked-in local fixture."""

    try:
        local_path = data_dir / "kahului_sample_1s.tt2"
        nc_path = tmp_path / "test.nc"

        # Write out NetCDF version of file
        ascii_topo = topotools.Topography(path=local_path)
        ascii_topo.read()
        ascii_topo.write(nc_path, topo_type=4, Z_format="%22.15e")

        # Read back in NetCDF file
        nc_topo = topotools.Topography(path=nc_path)
        nc_topo.read()

        assert np.allclose(ascii_topo.x, nc_topo.x), "Flat x-arrays did not match."
        assert np.allclose(ascii_topo.y, nc_topo.y), "Flat y-arrays did not match."
        assert np.allclose(ascii_topo.Z, nc_topo.Z), "Flat Z-arrays did not match."

    except ImportError:
        pytest.skip("Skipping test since NetCDF support not found.")
    except RuntimeError:
        pytest.skip("NetCDF topography test skipped due to runtime failure.")


# Remote/network integration test: keep opt-in and out of the default suite.
@pytest.mark.python
@pytest.mark.remote
def test_get_remote_file_remote(tmp_path):
    """Opt-in integration test for fetching a known remote topography file."""

    url = "".join((
        'https://raw.githubusercontent.com/rjleveque/geoclaw/',
        '5f675256c043e59e5065f9f3b5bdd41c2901702c/src/python/',
        'geoclaw/tests/kahului_sample_1s.tt2'
    ))
    try:
        clawpack.clawutil.data.get_remote_file(url, output_dir=tmp_path, force=True)
    except URLError:
        pytest.skip(f"Remote fetch failed for {url}. Skipping remote test.")
    
    local_path = tmp_path / Path(url).name
    assert local_path.exists(), (
        f"Expected file {local_path} not found after fetch attempt."
    )
    download_topo = topotools.Topography(path=local_path)

    test_path = data_dir / Path(url).name
    test_topo = topotools.Topography(path=test_path)

    assert download_topo.x.shape == test_topo.x.shape, (
        f"Downloaded x-array has wrong shape."
    )
    assert download_topo.y.shape == test_topo.y.shape,  (
        f"Downloaded y-array has wrong shape."
    )
    assert download_topo.Z.shape == test_topo.Z.shape,  (
        f"Downloaded Z-array has wrong shape."
    )
    assert np.allclose(download_topo.Z[:3, :3], test_topo.Z[:3, :3]), (
        f"Downloaded file does not match {test_path} in corner values."
    )
    assert np.allclose(download_topo.Z, test_topo.Z), (
        f"Downloaded file does not match {test_path}"
    )


# --- ETOPO1 integration tests and helpers ---

def _read_etopo1_topography(coarsen=10, return_xarray=False):
    """Read a small ETOPO1 subset for integration testing."""
    try:
        return topotools.read_netcdf(
            "etopo1",
            extent=etopo1_extent,
            coarsen=coarsen,
            return_xarray=return_xarray,
            verbose=True,
        )
    except (OSError, RuntimeError):
        pytest.skip("Reading ETOPO1 failed; check whether the remote server is available.")


@pytest.mark.python
@pytest.mark.netcdf
@pytest.mark.remote
def test_etopo1_topography():
    """Integration test for reading a remote ETOPO1 subset via topotools."""
    pytest.importorskip("netCDF4")

    topo1 = _read_etopo1_topography(coarsen=1)
    topo10 = _read_etopo1_topography(coarsen=10)

    testdata_path = data_dir / "etopo1_10min.asc"
    topo10input = topotools.Topography()
    topo10input.read(testdata_path, topo_type=3)

    assert topo1.Z.size > topo10.Z.size
    assert topo10.Z.shape == topo10input.Z.shape
    assert np.allclose(topo10.Z, topo10input.Z), (
        "topo10.Z does not agree with archived data"
    )


@pytest.mark.python
@pytest.mark.netcdf
@pytest.mark.remote
def test_etopo1_xarray():
    """Integration test for the xarray-returning ETOPO1 reader path."""
    pytest.importorskip("xarray")

    topo10, topo10_xarray = _read_etopo1_topography(coarsen=10, return_xarray=True)

    testdata_path = data_dir / "etopo1_10min.asc"
    topo10input = topotools.Topography()
    topo10input.read(testdata_path, topo_type=3)

    assert topo10.Z.shape == topo10input.Z.shape
    assert topo10_xarray["z"].shape == topo10input.Z.shape
    assert np.allclose(topo10_xarray["z"], topo10input.Z), (
        "topo10_xarray['z'] does not agree with archived data"
    )


def _import_pyplot():
    """Import pyplot using a non-interactive backend for test-safe plotting."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _make_unstructured_topo():
    """Construct a representative unstructured topography for interpolation tests."""
    # Create random test data
    f = lambda x, y: x * (1 - x) * np.cos(4 * np.pi * x) \
                                 * np.sin(4 * np.pi * y**2) ** 2

    fill_topo = topotools.Topography()
    fill_topo.x = np.linspace(0, 1, 100)
    fill_topo.y = np.linspace(0, 1, 200)
    fill_topo.Z = f(fill_topo.X, fill_topo.Y)

    points = np.loadtxt(data_dir / "unstructured_points.txt")
    values = f(points[:, 0], points[:, 1])

    # Create topography object
    topo = topotools.Topography(unstructured=True)
    topo.x = points[:, 0]
    topo.y = points[:, 1]
    topo.z = values

    return fill_topo, topo


@pytest.mark.python
def test_unstructured_topo():
    """Test interpolation from unstructured points onto a regular grid."""
    # Check to see if scipy is available and skip if not, since it's required
    # for interpolation.
    pytest.importorskip("scipy")

    fill_topo, topo = _make_unstructured_topo()
    topo.interp_unstructured(fill_topo, extent=[0, 1, 0, 1], delta=(1e-2, 1e-2))
    assert not topo.unstructured
    assert np.isfinite(topo.Z).all()
    
    compare_scalar = topo.Z[50, 50]
    assert np.isfinite(compare_scalar)

    test_data_path = data_dir / "unstructured_test_data.tt3"
    compare_data = topotools.Topography(path=test_data_path)
    assert np.allclose(compare_data.Z[50, 50], compare_scalar)
    assert np.allclose(compare_data.Z, topo.Z)


def save_unstructured_test_data(output_dir):
    """Utility function to save unstructured interpolation test data."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fill_topo, topo = _make_unstructured_topo()
    topo.interp_unstructured(fill_topo, extent=[0, 1, 0, 1], delta=(1e-2, 1e-2))
    test_data_path = output_dir / "unstructured_test_data.tt3"
    topo.write(test_data_path, Z_format="%22.15e")


def plot_unstructured_topo_baseline(output_dir):
    """Create optional diagnostic plots for the unstructured interpolation."""
    plt = _import_pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fill_topo, topo = _make_unstructured_topo()
    topo.interp_unstructured(fill_topo, extent=[0, 1, 0, 1], delta=(1e-2, 1e-2))

    fig = plt.figure(figsize=(16, 6))
    axes = fig.add_subplot(1, 3, 1)
    fill_topo.plot(axes=axes)
    axes.set_title("True Field")

    axes = fig.add_subplot(1, 3, 2)
    topo.plot(axes=axes, region_extent=[0, 1, 0, 1])
    axes.set_title("Unstructured Field")

    axes = fig.add_subplot(1, 3, 3)
    topo.plot(axes=axes)
    axes.set_title("Interpolated Field")

    fig.savefig(output_dir / "unstructured_interpolation.png")
    plt.close(fig)


def _make_bowl_hill_topography():
    """Construct a representative bowl-hill topography for plotting tests."""
    topo = topotools.Topography(topo_func=topo_bowl_hill)
    topo.x = np.linspace(-1.5, 2.5, 101)
    topo.y = np.linspace(-1.0, 2.0, 76)
    return topo


def plot_topo_bowl_hill(output_dir):
    """Create optional diagnostic plots for the bowl-hill topography."""
    plt = _import_pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    topo = _make_bowl_hill_topography()

    fig, ax = plt.subplots()
    topo.plot(axes=ax)
    fig.savefig(output_dir / "bowl_hill.png")
    plt.close(fig)

    topo_crop = topo.crop([0.5, 1.5, 0.0, 2.0])
    fig, ax = plt.subplots()
    topo_crop.plot(axes=ax)
    ax.set_title("Cropped topography")
    fig.savefig(output_dir / "bowl_hill_crop.png")
    plt.close(fig)


@pytest.mark.python
def test_plot_topo_bowl_hill():
    """Smoke test Topography.plot on the bowl-hill example."""
    plt = _import_pyplot()
    topo = _make_bowl_hill_topography()

    fig, ax = plt.subplots()
    topo.plot(axes=ax)
    assert ax.has_data()
    plt.close(fig)

    topo_crop = topo.crop([0.5, 1.5, 0.0, 2.0])
    fig, ax = plt.subplots()
    topo_crop.plot(axes=ax)
    assert ax.has_data()
    plt.close(fig)


def plot_kahului(output_dir):
    r"""Create optional diagnostic plots for the Kahului sample topography."""
    plt = _import_pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    path = data_dir / "kahului_sample_1s.tt2"
    topo = topotools.Topography(path, topo_type=2)

    fig, ax = plt.subplots()
    topo.plot(axes=ax)
    ax.set_title("Kahului Harbor at 1 second resolution")
    fig.savefig(output_dir / "kahului_imshow.png")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.contour(topo.X, topo.Y, topo.Z, np.linspace(-20, -2, 10), colors="b", linestyles="-")
    ax.contour(topo.X, topo.Y, topo.Z, np.linspace(2, 20, 10), colors="g")
    ax.contour(topo.X, topo.Y, topo.Z, [0.0], colors="r")

    mean_lat = 0.5 * (topo.y.max() + topo.y.min())
    ax.set_aspect(1.0 / np.cos(np.pi / 180.0 * mean_lat))
    ax.ticklabel_format(style="plain", useOffset=False)
    plt.xticks(rotation=20)
    ax.set_title("2-meter contours of topo (green) and bathymetry (blue)")
    fig.savefig(output_dir / "kahului_contour.png")
    plt.close(fig)


def plot_etopo1(output_dir):
    """Create optional diagnostic plots for the remote ETOPO1 integration tests."""
    plt = _import_pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pytest.importorskip("netCDF4")
    topo1 = _read_etopo1_topography(coarsen=1)
    topo10 = _read_etopo1_topography(coarsen=10)

    fig = plt.figure(figsize=(12, 5))
    ax1 = fig.add_subplot(1, 2, 1)
    topo1.plot(axes=ax1)
    ax1.set_title("1 minute ETOPO1 data")

    ax2 = fig.add_subplot(1, 2, 2)
    topo10.plot(axes=ax2)
    ax2.set_title("10 minute ETOPO1 data")

    fig.savefig(output_dir / "etopo1_test_plot.png")
    plt.close(fig)


@pytest.mark.python
def test_plot_kahului():
    r"""Smoke test plotting for a file-backed Topography object."""
    plt = _import_pyplot()

    path = data_dir / "kahului_sample_1s.tt2"
    topo = topotools.Topography(path, topo_type=2)

    assert topo.Z.shape == (46, 65), "*** K.Z is wrong shape"
    assert np.allclose(
        topo.Z[:3, :3],
        np.array(
            [
                [11.339, 11.339, 11.339],
                [13.339, 11.339, 11.339],
                [13.339, 11.339, 10.339],
            ]
        ),
    ), "*** Topography K does not match"

    fig, ax = plt.subplots()
    topo.plot(axes=ax)
    assert ax.has_data()
    plt.close(fig)

    fig, ax = plt.subplots()
    contour_sets = []
    contour_sets.append(ax.contour(topo.X, topo.Y, topo.Z, np.linspace(-20, -2, 10), colors="b", linestyles="-"))
    contour_sets.append(ax.contour(topo.X, topo.Y, topo.Z, np.linspace(2, 20, 10), colors="g"))
    contour_sets.append(ax.contour(topo.X, topo.Y, topo.Z, [0.0], colors="r"))
    assert all(len(cs.levels) > 0 for cs in contour_sets)

    mean_lat = 0.5 * (topo.y.max() + topo.y.min())
    aspect = 1.0 / np.cos(np.pi / 180.0 * mean_lat)
    assert np.isfinite(aspect) and aspect > 0.0
    ax.set_aspect(aspect)
    ax.ticklabel_format(style="plain", useOffset=False)
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == "save":
            output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".") / "unstructured_test_data"
            save_unstructured_test_data(output_dir)
        elif sys.argv[1].lower() == "plot":
            output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".") / "plot_output"
            plot_unstructured_topo_baseline(output_dir)
            plot_topo_bowl_hill(output_dir)
            plot_kahului(output_dir)
            plot_etopo1(output_dir)
        else:
            print("Usage: python test_topotools.py [save|plot] [output_dir]")
            print("Run remote tests via pytest, e.g.: pytest -m remote tests/test_topotools.py")
    else:
        raise SystemExit(pytest.main([__file__]))
