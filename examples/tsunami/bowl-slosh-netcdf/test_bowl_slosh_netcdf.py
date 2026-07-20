"""Pytest regression test for the GeoClaw bowl-slosh NetCDF example."""

from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

def _make_bowl_netcdf_topography(output_dir: Path, variant: str = "standard") -> None:
    """
    Create the NetCDF topography input used by the example.

    We are directly creating the NetCDF file here rather than relying on the
    example's maketopo.py to test out a slightly generic NetCDF writing approach
    that doesn't rely on the topotools write method.™

    Parameters
    ----------
    output_dir : Path
        Directory where the NetCDF file is written.
    variant : str
        Coordinate layout to produce.  ``"standard"`` (default) preserves the
        original behavior and writes ``bowl.nc``.  All other variants write a
        file named after the variant (see ``_VARIANT_FILENAMES``).
    """
    netCDF4 = pytest.importorskip("netCDF4")

    a = 1.0
    h0 = 0.1
    topo_func = lambda x, y: h0 * (x**2 + y**2) / a**2 - h0

    topo = topotools.Topography(topo_func=topo_func)
    topo.x = np.linspace(-3.1, 3.1, 310)
    topo.y = np.linspace(-3.5, 2.5, 300)

    filename = "bowl.nc"

    if variant == "standard":
        # Intentionally create latitude first to exercise dimension discovery.
        with netCDF4.Dataset(output_dir / filename, "w") as out:
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))
            elevations.units = "m"

            latitudes[:] = topo.y
            longitudes[:] = topo.x
            elevations[:] = topo.Z

    elif variant == "dim_order":
        # Data variable stored as (lon, lat) instead of (lat, lon); exercises
        # dimension-order discovery.
        with netCDF4.Dataset(output_dir / filename, "w") as out:
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lon", "lat"))
            elevations.units = "m"

            latitudes[:] = topo.y
            longitudes[:] = topo.x
            elevations[:] = topo.Z.T          # transpose to match (lon, lat)

    elif variant == "cf_compliant":
        # Standard (lat, lon) layout with full CF attributes added via
        # CFNormalizer; exercises CF-attribute-based coordinate discovery.
        from clawpack.geoclaw.netcdf_utils import CFNormalizer
        import xarray as xr

        ds = xr.Dataset(
            {"elevation": (["lat", "lon"], topo.Z, {"units": "m"})},
            coords={"lat": topo.y, "lon": topo.x},
        )
        ds_norm = CFNormalizer(ds).normalize()
        ds_norm.to_netcdf(
            output_dir / filename,
            encoding={"elevation": {"_FillValue": 9.96921e+36}},
        )

    elif variant == "cropped":
        # Topography domain is exactly 0.5 degrees larger than the simulation
        # extent on each side: x in [-2.5, 2.5], y in [-3.5, 2.5].
        # (Simulation domain after test setup: x in [-2, 2], y in [-3, 2].)
        topo_crop = topotools.Topography(topo_func=topo_func)
        topo_crop.x = np.linspace(-2.5, 2.5, 250)
        topo_crop.y = np.linspace(-3.5, 2.5, 300)

        with netCDF4.Dataset(output_dir / filename, "w") as out:
            out.createDimension("lat", len(topo_crop.y))
            out.createDimension("lon", len(topo_crop.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))
            elevations.units = "m"

            latitudes[:] = topo_crop.y
            longitudes[:] = topo_crop.x
            elevations[:] = topo_crop.Z

    elif variant == "coarsened":
        # Twice the standard resolution, corrupted at every point a coarsen=2
        # subsample should skip (odd row/column indices); the kept even-even
        # lattice reproduces the standard grid exactly.  Exercises the strided
        # NetCDF read path (tp_coarsen for topo_type=4).  The test pins the
        # kept window to index 0 with a full-extent crop_extent, since the AMR
        # domain fallback would start the window at an arbitrary parity.
        x = np.linspace(-3.1, 3.1, 619)
        y = np.linspace(-3.5, 2.5, 599)
        X, Y = np.meshgrid(x, y)
        Z = topo_func(X, Y)
        corrupt = np.ones(Z.shape, dtype=bool)
        corrupt[::2, ::2] = False
        Z[corrupt] -= 0.05

        with netCDF4.Dataset(output_dir / filename, "w") as out:
            out.createDimension("lat", len(y))
            out.createDimension("lon", len(x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))
            elevations.units = "m"

            latitudes[:] = y
            longitudes[:] = x
            elevations[:] = Z

    elif variant == "km":
        # Same bathymetry, but elevation stored in kilometers (units='km').
        # The descriptor scale_factor (1000) must convert it to meters on the
        # Fortran read, reproducing the meters baseline exactly.  Exercises the
        # Fortran-side unit scaling for topo_type=4.
        with netCDF4.Dataset(output_dir / filename, "w") as out:
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))
            elevations.units = "km"

            latitudes[:] = topo.y
            longitudes[:] = topo.x
            elevations[:] = topo.Z / 1000.0

    else:
        raise ValueError(f"Unknown variant '{variant}'.")


@pytest.mark.regression
@pytest.mark.netcdf
@pytest.mark.parametrize("variant", ["standard", "dim_order", "cf_compliant",
                                     "cropped", "coarsened", "km"])
def test_bowl_slosh_netcdf(tmp_path: Path, save: bool, variant: str) -> None:
    """
    Parametrized regression test exercising non-standard NetCDF coordinate layouts.

    Each variant writes a topography file with a different coordinate convention
    while encoding the same parabolic-bowl bathymetry as test_bowl_slosh_netcdf.
    Because the physical bathymetry is identical, all variants must produce gauge
    output that matches the existing reference data.

    Variants
    --------
    standard     Intentionally create latitude first to exercise dimension
                 discovery; otherwise, a pretty standard layout.
    dim_order    Data variable stored as (lon, lat) instead of (lat, lon);
                 exercises dimension-order discovery.
    cf_compliant Standard (lat, lon) layout with full CF attributes added via
                 CFNormalizer (standard_name, axis, units, _FillValue); exercises
                 CF-attribute-based coordinate discovery.
    cropped      Topography domain is exactly 0.5 degrees larger than the
                 simulation extent on each side (x in [-2.5, 2.5],
                 y in [-3.5, 2.5]); exercises minimal-margin domain coverage.
    coarsened    Twice the standard resolution with corrupted values at the
                 points a coarsen=2 subsample must skip; exercises the strided
                 NetCDF read (tp_coarsen for topo_type=4) with detection power
                 against a missed or mis-indexed stride.
    """
    pytest.importorskip("netCDF4")

    _make_bowl_netcdf_topography(tmp_path, variant=variant)

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    runner.set_data()
    runner.rundata.clawdata.lower[1] = -3.0
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.5
    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append([1, 0.5, 0.5, 0.0, 1e10])

    if variant == "coarsened":
        topo_entry = topotools.Topography()
        topo_entry.path = "bowl.nc"
        topo_entry.topo_type = 4
        topo_entry.coarsen = 2
        # Full-extent crop pins the kept lattice to file index 0 (see the
        # variant writer); the subsample then lands on the standard grid.
        topo_entry.crop_extent = [-3.1, 3.1, -3.5, 2.5]
        runner.rundata.topo_data.topofiles = [topo_entry]

    runner.write_data()

    runner.build_executable()
    runner.run_code()

    runner.check_gauge(gauge_id=1, indices=(2, 3), rtol=1.0e-4, atol=1.0e-4, save=save)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
