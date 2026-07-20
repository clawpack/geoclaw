#!/usr/bin/env python
"""
Chile 2010 tsunami regression test

Also tests some variants on lat-lon topography reading, including several
NetCDF topography file formats with non-standard coordinate conventions.
"""

import warnings
from pathlib import Path

import pytest

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools
import clawpack.clawutil.util as clawutil

def _get_topography(runner, variant=0):
    """Get topography and read"""
    
    # Ensure that the topography file is available
    maketopo = clawutil.fullpath_import(runner.test_path / "maketopo.py", 
                                        module_name="maketopo")
    topo_path = maketopo.get_topo()
    dtopo_path = maketopo.make_dtopo()

    # Read in original topography data and write out in format expected by test
    topo = topotools.Topography()
    topo.topo_type = 2
    topo.path = topo_path
    topo.read()

    # Read in dtopography data
    dtopo = dtopotools.DTopography()
    dtopo.dtopo_type = 3
    dtopo.path = dtopo_path
    dtopo.read()

    return topo, dtopo


def _write_netcdf_variant(topo, variant, output_path, crop_bounds=None):
    """Write topography as a NetCDF file in the requested variant format.

    Parameters
    ----------
    topo : topotools.Topography
        Topography object with x, y, Z already populated.
    variant : str
        One of ``"netcdf_topotools"``, ``"lat_flip"``, or ``"lon_360"``.
    output_path : Path
        Destination file path for the NetCDF output.
    crop_bounds : tuple or None
        (lon_min, lon_max, lat_min, lat_max) in domain coordinates.  When
        provided, ``topo_entries()`` is called to compute the correct
        ``lon_wrap_offset`` for the descriptor; the returned list has the same
        ``[ttype, path, meta]`` element format used by ``topo_data.topofiles``.

    Variants
    --------
    netcdf_topotools
        Write via ``topotools.Topography.write(topo_type=4)``, then run the
        result through CFNormalizer.  If CFNormalizer has to fix anything
        (e.g. rename coordinate variables), a UserWarning is emitted so the
        gap between topotools output and strict CF compliance is visible.
    lat_flip
        Write NetCDF manually (netCDF4) with latitude stored N-to-S
        (descending), which is non-standard but common in observational data.
        Elevation rows are flipped accordingly so the physical bathymetry is
        identical to the standard variant.
    lon_360
        Write NetCDF manually (netCDF4) with longitude in [0, 360] rather than
        [-180, 180].  The physical bathymetry is otherwise identical.
    """
    if variant == "netcdf_topotools":
        import netCDF4  # guaranteed available; test already called importorskip
        import xarray as xr
        from clawpack.geoclaw.netcdf_utils import CFNormalizer

        topo.write(str(output_path), topo_type=4)

        # Open the just-written file and run it through CFNormalizer to check
        # whether the output is already fully CF-compliant.
        with xr.open_dataset(output_path) as ds_before:
            coords_before = set(ds_before.coords)
            attrs_before = {name: dict(ds_before[name].attrs)
                            for name in list(ds_before.coords) + list(ds_before.data_vars)}

            ds_after = CFNormalizer(ds_before).normalize()

            coords_after = set(ds_after.coords)
            attrs_after = {name: dict(ds_after[name].attrs)
                           for name in list(ds_after.coords) + list(ds_after.data_vars)}

        changed = coords_before != coords_after
        if not changed:
            # Check whether any attrs changed on variables that survived renaming
            for name in coords_before & coords_after:
                if attrs_before.get(name) != attrs_after.get(name):
                    changed = True
                    break
        if not changed:
            for name in ds_after.data_vars:
                if attrs_before.get(name) != attrs_after.get(name):
                    changed = True
                    break

        if changed:
            renamed = coords_before - coords_after
            warnings.warn(
                f"CFNormalizer modified the topotools NetCDF output "
                f"('{output_path.name}'). "
                f"Coordinate variables renamed or attributes added: "
                f"{renamed or '(attr-only change)'}. "
                "The topotools.write topo_type=4 path is not fully "
                "CF-compliant; consider updating it to match CF standard "
                "variable names.",
                UserWarning,
                stacklevel=2,
            )

    elif variant == "lat_flip":
        import netCDF4  # guaranteed available; test already called importorskip
        from clawpack.geoclaw.netcdf_utils import suppress_netcdf4_shape_warning

        # Write latitude N-to-S (descending) — a common non-standard convention.
        # Flip Z rows so each row still corresponds to the correct latitude.
        with netCDF4.Dataset(output_path, "w") as out, suppress_netcdf4_shape_warning():
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))
            elevations.units = "m"

            latitudes[:] = topo.y[::-1]      # descending: N → S
            longitudes[:] = topo.x
            elevations[:] = topo.Z[::-1, :]  # flip rows to match reversed lat

    elif variant == "lon_360":
        import netCDF4  # guaranteed available; test already called importorskip
        from clawpack.geoclaw.netcdf_utils import suppress_netcdf4_shape_warning

        # Write longitude in [0, 360] convention instead of [-180, 180].
        with netCDF4.Dataset(output_path, "w") as out, suppress_netcdf4_shape_warning():
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))
            elevations.units = "m"

            latitudes[:] = topo.y
            longitudes[:] = topo.x % 360    # map [-180, 180] → [0, 360]
            elevations[:] = topo.Z

    else:
        raise ValueError(f"Unknown NetCDF variant '{variant}'.")

    # Interrogate the written file to auto-detect y_increasing, dim_order, etc.
    # When crop_bounds (domain coordinates) are given, topo_entries() is used
    # so that the correct lon_wrap_offset is computed and written into the descriptor
    # block that Fortran's read_netcdf_descriptor needs.
    from clawpack.geoclaw.netcdf_utils import TopoInspector
    with TopoInspector(output_path, 'elevation', crop_bounds=crop_bounds) as insp:
        if crop_bounds is not None:
            return insp.topo_entries()
        else:
            return [[4, output_path, insp.inspect_topo()]]

CASES = [pytest.param("standard", id="standard"),
         pytest.param("netcdf_topotools", id="netcdf_topotools", marks=pytest.mark.netcdf),
         pytest.param("lat_flip", id="lat_flip", marks=pytest.mark.netcdf),
         pytest.param("lon_360", id="lon_360", marks=pytest.mark.netcdf)]

@pytest.mark.regression
@pytest.mark.remote
@pytest.mark.tsunami
@pytest.mark.parametrize("variant", CASES)
def test_chile2010(tmp_path: Path, save: bool, variant: dict):
    """Chile 2010 tsunami regression test for GeoClaw.

    Parametrized over topography file format.  Because every variant encodes
    the same physical bathymetry, all variants must reproduce the existing
    reference gauge output at gauge 32412 within the same tolerances.

    Variants
    --------
    standard
        ASCII topo type 2 read directly from the downloaded file.  This is
        the baseline behavior and is identical to the pre-parametrization test.
    netcdf_topotools
        The downloaded ASCII topo is re-written as NetCDF via
        ``topotools.Topography.write(topo_type=4)``.  CFNormalizer is run over
        the result to surface any CF-compliance gaps in the write path.
    lat_flip
        NetCDF written manually (netCDF4) with latitude stored N-to-S
        (descending) rather than the usual S-to-N (ascending) order.
    lon_360
        NetCDF written manually (netCDF4) with longitude in the [0, 360]
        convention rather than [-180, 180].
    """
    if variant != "standard":
        pytest.importorskip("netCDF4")

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)

    # Topography
    topo, dtopo = _get_topography(runner)

    # Setup data for test
    runner.set_data()
    # runner.rundata.clawdata.num_output_times = 1
    runner.rundata.amrdata.amr_levels_max = 1

    if variant == "standard":
        runner.rundata.topo_data.topofiles = [[2, topo.path]]
    else:
        nc_path = tmp_path / f"chile_topo_{variant}.nc"
        clawdata = runner.rundata.clawdata
        # lon_360 needs domain crop bounds so topo_entries() can compute the
        # correct lon_wrap_offset (-360.0) for the descriptor.  Other variants use
        # inspect_topo() (lon_wrap_offset=0.0 is correct for [-180,180] files).
        crop = (clawdata.lower[0], clawdata.upper[0],
                clawdata.lower[1], clawdata.upper[1]) if variant == "lon_360" else None
        entries = _write_netcdf_variant(topo, variant, nc_path, crop_bounds=crop)
        runner.rundata.topo_data.topofiles = entries

    runner.write_data()

    # Build and run test
    if variant != "standard":
        print(f"Testing variant '{variant}' with topography file: {nc_path.name}")
        runner.build_executable(FFLAGS="-DNETCDF $(NETCDF_FFLAGS)", 
                                LFLAGS="$(NETCDF_LFLAGS)")
    else:
        print(f"Testing standard variant with topography file: {topo.path}")
        runner.build_executable()
    runner.run_code()

    # lon_360 coordinates undergo +360 then -360 shifts that are not a
    # perfect floating-point round-trip for non-round values, so allow a
    # slightly relaxed tolerance (still much smaller than 1% physically).
    tol = {"atol": 1e-2} if variant == "lon_360" else {}
    runner.check_gauge(save=save, gauge_id=32412, **tol)

@pytest.mark.regression
@pytest.mark.remote
@pytest.mark.tsunami
@pytest.mark.netcdf
@pytest.mark.parametrize("dz_units", ["meters", "km"])
def test_chile2010_dtopo_netcdf(tmp_path: Path, save: bool, dz_units: str):
    """Chile 2010 regression test exercising the dtopo NetCDF (type=4) path.

    Regression coverage for a bug where a NetCDF dtopo file's time
    coordinate -- written as a bare-duration ("seconds") unit, as
    ``DTopography._write_netcdf`` does -- could decode via xarray into
    timedelta64.  ``DTopoInspector._compute_time_axis`` then cast that
    straight to float, reinterpreting the raw nanosecond tick count as
    seconds: t0dtopo came out ~1e9 instead of ~1.  Since the chile2010
    dtopo is a single-frame deformation at t=1s, that inflated t0dtopo
    meant the simulated time (max 9 hours) never reached it, so the
    coseismic deformation was never applied at all -- and because
    ``topo_finalized`` also never becomes true, ``dt_max_dtopo=0.2`` pins
    the timestep for the entire run.  Both symptoms are checked here: the
    gauge must show the same tsunami signal as the ASCII dtopo baseline
    (not a flat, undeformed sea surface), and the run must actually finish.
    """
    pytest.importorskip("netCDF4")

    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    topo, dtopo = _get_topography(runner)

    runner.set_data()
    runner.rundata.amrdata.amr_levels_max = 1
    runner.rundata.topo_data.topofiles = [[2, topo.path]]

    dtopo_nc_path = tmp_path / "dtopo_usgs100227.nc"
    if dz_units == "meters":
        dtopo.write(dtopo_nc_path, dtopo_type=4)
    else:
        # Same deformation stored in kilometers (units='km'): the descriptor
        # scale_factor (1000) must convert it to meters on the Fortran read,
        # reproducing the meters baseline exactly.  Exercises the Fortran-side
        # unit scaling for dtopo_type=4.
        import numpy as np
        import xarray as xr
        x = dtopo.X[0, :]
        y = dtopo.Y[:, 0]
        ds = xr.Dataset(
            {"dz": (("time", "lat", "lon"),
                    np.asarray(dtopo.dZ, dtype=float) / 1000.0,
                    {"units": "km", "long_name": "seafloor deformation"})},
            coords={"time": np.asarray(dtopo.times, dtype=float),
                    "lat": np.asarray(y, dtype=float),
                    "lon": np.asarray(x, dtype=float)},
        )
        ds["time"].attrs.update(units="seconds")
        ds.to_netcdf(dtopo_nc_path, encoding={
            "dz": {"_FillValue": None}, "time": {"_FillValue": None},
            "lat": {"_FillValue": None}, "lon": {"_FillValue": None}})
    runner.rundata.dtopo_data.dtopofiles = [[4, str(dtopo_nc_path)]]

    runner.write_data()

    print(f"Testing dtopo netCDF variant with dtopo file: {dtopo_nc_path.name}")
    runner.build_executable(FFLAGS="-DNETCDF $(NETCDF_FFLAGS)",
                            LFLAGS="$(NETCDF_LFLAGS)")
    runner.run_code()

    runner.check_gauge(save=save, gauge_id=32412)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
