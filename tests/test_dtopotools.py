#!/usr/bin/env python

import sys
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.dtopotools as dtopotools

# Local test directory and bundled test data
testdir = Path(__file__).parent
data_dir = testdir / "data"


def _make_grid(xlower, xupper, ylower, yupper, points_per_degree=4):
    """Construct a regular longitude-latitude grid for dtopography tests."""
    dx = 1.0 / points_per_degree
    mx = int((xupper - xlower) / dx + 1)
    xupper = xlower + (mx - 1) * dx
    my = int((yupper - ylower) / dx + 1)
    yupper = ylower + (my - 1) * dx
    x = np.linspace(xlower, xupper, mx)
    y = np.linspace(ylower, yupper, my)
    return x, y


def _assert_dtopo_matches_baseline(dtopo, test_data_path):
    """Compare a generated DTopography object against a bundled baseline."""
    compare_data = dtopotools.DTopography(path=test_data_path)
    compare_data.read(path=test_data_path, dtopo_type=3)

    assert dtopo.dZ.shape == compare_data.dZ.shape, (
        f"dtopo.dZ.shape is {dtopo.dZ.shape}, should be {compare_data.dZ.shape}"
    )
    assert np.allclose(compare_data.dZ, dtopo.dZ)


def _import_pyplot():
    """Import pyplot using a non-interactive backend for test-safe plotting."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def save_dtopo_test_data(output_dir):
    """Utility helper to regenerate bundled dtopography baselines."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_units = {
        "length": "km",
        "width": "km",
        "depth": "km",
        "slip": "m",
        "mu": "dyne/cm^2",
    }

    csv_fault = dtopotools.CSVFault()
    csv_fault.read(
        data_dir / "alaska1964.csv",
        input_units=input_units,
        coordinate_specification="noaa sift",
    )
    x, y = _make_grid(203.0, 214.0, 55.0, 60.0)
    csv_fault.create_dtopography(x, y, times=[1.0]).write(
        output_dir / "alaska1964_test_data.tt3", dtopo_type=3
    )

    ucsb_fault = dtopotools.UCSBFault()
    ucsb_fault.read(data_dir / "tohoku_ucsb.txt")
    x, y = _make_grid(140.0, 146.0, 35.0, 41.0)
    tmax = 0.0
    for subfault in ucsb_fault.subfaults:
        tmax = max(tmax, subfault.rupture_time + subfault.rise_time)
    ucsb_fault.rupture_type = "kinematic"
    times = np.linspace(0.0, tmax, 10)
    ucsb_fault.create_dtopography(x, y, times).write(
        output_dir / "tohoku_test_data.tt3", dtopo_type=3
    )

    sift_fault = dtopotools.SiftFault({"acsza1": 2, "acszb1": 3})
    x, y = _make_grid(162.0, 168.0, 53.0, 59.0)
    sift_fault.create_dtopography(x, y, [1.0]).write(
        output_dir / "sift_test_data.tt3", dtopo_type=3
    )

    subdivided_fault = dtopotools.SiftFault({"acsza1": 1.0})
    fault_plane = subdivided_fault.subfaults[0]
    subdivided = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    x, y = _make_grid(162.0, 168.0, 53.0, 59.0)
    subdivided.create_dtopography(x, y, [1.0]).write(
        output_dir / "SubdividedFaultPlane_test_data.tt3", dtopo_type=3
    )


@pytest.mark.python
def test_read_csv_make_dtopo():
    r"""Test reading and making of a CSV subfault specified dtopo."""
    subfault_path = data_dir / "alaska1964.csv"
    input_units = {
        "length": "km",
        "width": "km",
        "depth": "km",
        "slip": "m",
        "mu": "dyne/cm^2",
    }
    fault = dtopotools.CSVFault()
    fault.read(
        subfault_path,
        input_units=input_units,
        coordinate_specification="noaa sift",
    )

    assert abs(fault.Mw() - 8.53336) < 1e-4, f"*** Mw is wrong: {fault.Mw():g}"

    x, y = _make_grid(203.0, 214.0, 55.0, 60.0)
    dtopo = fault.create_dtopography(x, y, times=[1.0])

    test_data_path = data_dir / "alaska1964_test_data.tt3"
    _assert_dtopo_matches_baseline(dtopo, test_data_path)


@pytest.mark.python
def test_read_ucsb_make_dtopo():
    r"""Test reading and making of a UCSB subfault specified dtopo."""
    subfault_path = data_dir / "tohoku_ucsb.txt"
    fault = dtopotools.UCSBFault()
    fault.read(subfault_path)

    assert abs(fault.Mw() - 9.13957973) < 1e-4, f"*** Mw is wrong: {fault.Mw():g}"

    x, y = _make_grid(140.0, 146.0, 35.0, 41.0)

    tmax = 0.0
    for subfault in fault.subfaults:
        tmax = max(tmax, subfault.rupture_time + subfault.rise_time)

    fault.rupture_type = "kinematic"
    times = np.linspace(0.0, tmax, 10)
    dtopo = fault.create_dtopography(x, y, times)

    test_data_path = data_dir / "tohoku_test_data.tt3"
    _assert_dtopo_matches_baseline(dtopo, test_data_path)


@pytest.mark.python
def test_read_sift_make_dtopo():
    r"""Test reading and making of a SIFT subfault specified dtopo."""
    sift_slip = {"acsza1": 2, "acszb1": 3}
    fault = dtopotools.SiftFault(sift_slip)

    assert abs(fault.Mw() - 7.966666666667) < 1e-4, f"*** Mw is wrong: {fault.Mw():g}"

    x, y = _make_grid(162.0, 168.0, 53.0, 59.0)
    dtopo = fault.create_dtopography(x, y, [1.0])

    test_data_path = data_dir / "sift_test_data.tt3"
    _assert_dtopo_matches_baseline(dtopo, test_data_path)


@pytest.mark.python
def test_subdivided_plane_fault_make_dtopo():
    r"""Test dtopography generation from a subdivided plane fault."""
    sift_slip = {"acsza1": 1.0}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]

    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)

    assert abs(fault2.Mw() - 7.500686667) < 1e-4, (
        f"*** Mw is wrong: {fault2.Mw():g}"
    )

    x, y = _make_grid(162.0, 168.0, 53.0, 59.0)
    dtopo = fault2.create_dtopography(x, y, [1.0])

    test_data_path = data_dir / "SubdividedFaultPlane_test_data.tt3"
    _assert_dtopo_matches_baseline(dtopo, test_data_path)


@pytest.mark.python
@pytest.mark.parametrize(
    "filename, dtopo_type",
    [
        ("alaska1964.tt1", 1),
        ("alaska1964.tt3", 3),
    ],
)
def test_dtopo_io(tmp_path, filename, dtopo_type):
    # Written with full floating-point precision, the write/read round-trip is
    # exact; any loss seen with the default ``dZ_format="%.3f"`` is text-format
    # rounding (~5e-4), not a defect in the reader/writer.
    test_data_path = data_dir / "alaska1964_test_data.tt3"
    test_dtopo = dtopotools.DTopography(path=test_data_path)
    test_dtopo.read(path=test_data_path, dtopo_type=3)

    path = tmp_path / filename
    test_dtopo.write(path, dtopo_type=dtopo_type, dZ_format="%.12e")

    dtopo = dtopotools.DTopography()
    dtopo.read(path=path, dtopo_type=dtopo_type)

    assert test_dtopo.dZ.shape == dtopo.dZ.shape
    assert np.allclose(test_dtopo.dZ, dtopo.dZ, atol=1e-10, rtol=0.0)


def _make_synthetic_dtopo():
    """Construct a small two-frame DTopography with a Gaussian uplift."""
    x, y = _make_grid(0.0, 1.0, 0.0, 1.0)
    X, Y = np.meshgrid(x, y)
    dtopo = dtopotools.DTopography()
    dtopo.x = x
    dtopo.y = y
    dtopo.X = X
    dtopo.Y = Y
    dtopo.times = [0.0, 1.0]
    dZ_final = np.exp(-((X - 0.5)**2 + (Y - 0.5)**2) / 0.1)
    dtopo.dZ = np.array([np.zeros(X.shape), dZ_final])
    return dtopo


@pytest.mark.python
def test_dtopo_preprocessing_attributes(tmp_path):
    r"""negate_z and z_shift are applied by DTopography.read()."""
    path = tmp_path / "synthetic.tt3"
    _make_synthetic_dtopo().write(path, dtopo_type=3, dZ_format="%.12e")

    plain = dtopotools.DTopography()
    plain.read(path=path, dtopo_type=3)

    negated = dtopotools.DTopography()
    negated.negate_z = True
    negated.read(path=path, dtopo_type=3)
    assert np.allclose(negated.dZ, -plain.dZ)

    shifted = dtopotools.DTopography()
    shifted.z_shift = 0.5
    shifted.read(path=path, dtopo_type=3)
    assert np.allclose(shifted.dZ, plain.dZ + 0.5)

    # Order matches topotools: negate first, then shift.
    both = dtopotools.DTopography()
    both.negate_z = True
    both.z_shift = 0.5
    both.read(path=path, dtopo_type=3)
    assert np.allclose(both.dZ, -plain.dZ + 0.5)


@pytest.mark.python
@pytest.mark.netcdf
def test_dtopo_netcdf_roundtrip(tmp_path):
    r"""dtopo_type=4 write/read round-trips dZ, times, and coordinates."""
    pytest.importorskip("netCDF4")

    base = _make_synthetic_dtopo()
    path = tmp_path / "synthetic.nc"
    base.write(path, dtopo_type=4)

    dtopo = dtopotools.DTopography()
    dtopo.read(path=path)   # type 4 inferred from .nc extension

    assert dtopo.dZ.shape == base.dZ.shape
    assert np.allclose(dtopo.dZ, base.dZ)
    assert np.allclose(dtopo.times, base.times)
    assert np.allclose(dtopo.x, base.x)
    assert np.allclose(dtopo.y, base.y)

    # Preprocessing attributes also apply on the NetCDF read path.
    shifted = dtopotools.DTopography()
    shifted.negate_z = True
    shifted.z_shift = 0.5
    shifted.read(path=path)
    assert np.allclose(shifted.dZ, -base.dZ + 0.5)


@pytest.mark.python
@pytest.mark.parametrize("attr, value", [
    ("crop_extent", [0.0, 1.0, 0.0, 1.0]),
    ("coarsen", 2),
    ("buffer", 1.0),
    ("align", (0.0, 0.5)),
    ("x_shift", 0.5),
])
def test_dtopo_unsupported_preprocessing(tmp_path, attr, value):
    r"""Unsupported preprocessing attributes raise in DTopography.read()."""
    path = tmp_path / "synthetic.tt3"
    _make_synthetic_dtopo().write(path, dtopo_type=3, dZ_format="%.12e")

    dtopo = dtopotools.DTopography()
    setattr(dtopo, attr, value)
    with pytest.raises(NotImplementedError, match=attr):
        dtopo.read(path=path, dtopo_type=3)


@pytest.mark.python
def test_geometry():
    r"""Test subfault geometry calculation."""

    pytest.skip(
        "Skipping geometry test for now; revisit whether it adds unique coverage beyond dtopography-generation tests."
    )
    from . import old_dtopotools

    subfault_path = data_dir / 'alaska1964.csv'
    input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m",
         "mu":"dyne/cm^2"}
    specifications = ['top center', 'centroid', 'bottom center', 'noaa sift']

    for specification in specifications:
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units,
                                  coordinate_specification=specification)

        # Subfault 10 is chosen at random, maybe do all?
        subfault = fault.subfaults[10]
        geometry = old_dtopotools.set_geometry(subfault)

        coord_tests = {
            "top center": {'test': [geometry['x_top'],
                                    geometry['y_top'],
                                    geometry['depth_top']],
                           'computed': subfault.centers[0]},
            "centroid": {'test': [geometry['x_centroid'],
                                  geometry['y_centroid']],
                         'computed': subfault.centers[1][:2]},
            "bottom center": {"test": [geometry['x_bottom'],
                                       geometry['y_bottom'],
                                       geometry['depth_bottom']],
                              "computed": subfault.centers[2]},
            "Corner A": {"test": [geometry["x_corners"][2],
                                  geometry["y_corners"][2]],
                         "computed": subfault.corners[0][:2]},
            "Corner B": {"test": [geometry["x_corners"][3],
                                  geometry["y_corners"][3]],
                         "computed": subfault.corners[1][:2]},
            "Corner C": {"test": [geometry["x_corners"][0],
                                  geometry["y_corners"][0]],
                         "computed": subfault.corners[2][:2]},
            "Corner D": {"test": [geometry["x_corners"][1],
                                  geometry["y_corners"][1]],
                         "computed": subfault.corners[3][:2]}
        }

        for (values, coord_test) in coord_tests.items():
            assert np.allclose(coord_test['test'], coord_test['computed']), (
                "Specification = %s, coords= %s:\n%s !=\n%s"
                % (specification, values, coord_test['test'], coord_test['computed'])
            )


@pytest.mark.python
def test_vs_old_dtopo():
    r"""Test new dtopotools with old version from 5.2."""
    pytest.skip("Skipping comparison with old tools.")


def _make_dynamic_tohoku_dtopo(verbose=False):
    """Construct the dynamic Tohoku dtopography test case."""
    shoreline_fname = data_dir / "tohoku_shoreline_1min.npy"
    shoreline_xy = np.load(shoreline_fname)

    subfault_fname = data_dir / "tohoku_ucsb.txt"
    fault = dtopotools.UCSBFault()
    fault.read(subfault_fname)
    fault.rupture_type = "dynamic"

    xlower = 140.0
    xupper = 146.0
    ylower = 35.0
    yupper = 41.0
    xylim = [xlower, xupper, ylower, yupper]

    mx = int((xupper - xlower) * 15 + 1)
    my = int((yupper - ylower) * 15 + 1)

    x = np.linspace(xlower, xupper, mx)
    y = np.linspace(ylower, yupper, my)

    tmax = 0.0
    for subfault in fault.subfaults:
        tmax = max(tmax, subfault.rupture_time + subfault.rise_time + subfault.rise_time_ending)
    if verbose:
        print("rupture ends at time", tmax)

    times = np.linspace(0.0, tmax, 10)
    dtopo = fault.create_dtopography(x, y, times, verbose=True)
    return shoreline_xy, fault, dtopo, xylim, tmax


@pytest.mark.python
def test_dynamic_tohoku(tmp_path):
    r"""Test dynamic faulting via a Tohoku example."""
    _, _, dtopo, _, _ = _make_dynamic_tohoku_dtopo()

    fname = tmp_path / "tohoku_ucsb_dynamic.tt3"
    dtopo.write(fname, 3)
    assert fname.exists()


def plot_dynamic_tohoku(output_dir):
    """Create optional diagnostic plots for the dynamic Tohoku example."""
    plt = _import_pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    shoreline_xy, fault, dtopo, xylim, tmax = _make_dynamic_tohoku_dtopo(verbose=True)
    dz_max = dtopo.dZ[-1].max()

    fig = plt.figure(figsize=(12, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fault.plot_subfaults(axes=ax1, slip_color=True, slip_time=tmax, xylim=xylim)
    dtopo.plot_dz_colors(axes=ax2, t=tmax, cmax_dz=dz_max)
    ax1.plot(shoreline_xy[:, 0], shoreline_xy[:, 1], "g")
    ax2.plot(shoreline_xy[:, 0], shoreline_xy[:, 1], "g")
    ax1.axis(xylim)
    ax2.axis(xylim)
    fig.savefig(output_dir / "tohoku_dynamic.png")
    plt.close(fig)


@pytest.mark.python
def test_subdivided_plane_fault():
    r"""Test the SubdividedPlaneFault class."""
    sift_slip = {"acsza1": 1.0}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]
    original_mo = fault_plane.Mo()

    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    assert fault2.Mo() > 0.0

    slip_function = lambda xi, eta: xi * (1 - xi) * eta
    fault2 = dtopotools.SubdividedPlaneFault(
        fault_plane,
        nstrike=5,
        ndip=3,
        slip_function=slip_function,
        Mo=original_mo,
    )
    assert np.isclose(fault2.Mo(), original_mo)

    fault2.subdivide(nstrike=20, ndip=10, slip_function=slip_function, Mo=original_mo)
    assert np.isclose(fault2.Mo(), original_mo)


def plot_subdivided_plane_fault(output_dir):
    """Create optional diagnostic plots for subdivided plane-fault examples."""
    plt = _import_pyplot()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sift_slip = {"acsza1": 1.0}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]
    original_mo = fault_plane.Mo()

    fig = plt.figure()
    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    fault2.plot_subfaults(slip_color=True)
    fig.savefig(output_dir / "subdivided_plane_fault_uniform.png")
    plt.close(fig)

    slip_function = lambda xi, eta: xi * (1 - xi) * eta
    fig = plt.figure()
    fault3 = dtopotools.SubdividedPlaneFault(
        fault_plane,
        nstrike=5,
        ndip=3,
        slip_function=slip_function,
        Mo=original_mo,
    )
    fault3.plot_subfaults(slip_color=True)
    fig.savefig(output_dir / "subdivided_plane_fault_weighted.png")
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "plot":
            output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".") / "plot_output"
            plot_dynamic_tohoku(output_dir)
            plot_subdivided_plane_fault(output_dir)
        elif sys.argv[1] == "save":
            output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else data_dir
            save_dtopo_test_data(output_dir)
        else:
            print("Usage: python test_dtopotools.py [save|plot] [output_dir]")
    else:
        raise SystemExit(pytest.main([__file__]))
