#!/usr/bin/env python
"""
Chile 2010 tsunami regression test

Also tests some variants on lat-lon topography reading, including several
NetCDF topography file formats with non-standard coordinate conventions.
"""

import warnings
from pathlib import Path

import pytest

import numpy as np

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools
import clawpack.clawutil.util as clawutil

def _get_topography(runner, variant=0):
    """Get topography and read"""

    # Ensure that the topography file is available
    maketopo = clawutil.fullpath_import("maketopo.py", module_name="maketopo")
    topo_path = maketopo.get_topo()

    # Read in original topography data and write out in format expected by test
    topo = topotools.Topography()
    topo.topo_type = 2
    topo.path = topo_path
    topo.read()

    return topo


def _write_netcdf_variant(topo, variant, output_path):
    """Write topography as a NetCDF file in the requested variant format.

    Parameters
    ----------
    topo : topotools.Topography
        Topography object with x, y, Z already populated.
    variant : str
        One of ``"netcdf_topotools"``, ``"lat_flip"``, or ``"lon_360"``.
    output_path : Path
        Destination file path for the NetCDF output.

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

        # Write latitude N-to-S (descending) — a common non-standard convention.
        # Flip Z rows so each row still corresponds to the correct latitude.
        with netCDF4.Dataset(output_path, "w") as out:
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))

            latitudes[:] = topo.y[::-1]      # descending: N → S
            longitudes[:] = topo.x
            elevations[:] = topo.Z[::-1, :]  # flip rows to match reversed lat

    elif variant == "lon_360":
        import netCDF4  # guaranteed available; test already called importorskip

        # Write longitude in [0, 360] convention instead of [-180, 180].
        with netCDF4.Dataset(output_path, "w") as out:
            out.createDimension("lat", len(topo.y))
            out.createDimension("lon", len(topo.x))

            latitudes = out.createVariable("lat", "f8", ("lat",))
            longitudes = out.createVariable("lon", "f8", ("lon",))
            elevations = out.createVariable("elevation", "f8", ("lat", "lon"))

            latitudes[:] = topo.y
            longitudes[:] = topo.x % 360    # map [-180, 180] → [0, 360]
            elevations[:] = topo.Z

    else:
        raise ValueError(f"Unknown NetCDF variant '{variant}'.")

    # Interrogate the written file to auto-detect lon_convention, lat_order,
    # and dim_order.  The returned TopoMetadata is passed to topofiles as the
    # third element so data.py writes the descriptor block that Fortran's
    # read_netcdf_descriptor needs.  Without it nc_lon_convention defaults to
    # 180 and [0,360] longitudes are never shifted back to [-180,180], causing
    # the "topo arrays do not cover domain" check to fail.
    from clawpack.geoclaw.netcdf_utils import TopoInterrogator
    return TopoInterrogator(output_path, 'elevation').interrogate_topo()


@pytest.mark.regression
@pytest.mark.remote
@pytest.mark.tsunami
@pytest.mark.parametrize("variant", ["standard", "netcdf_topotools", "lat_flip", "lon_360"])
def test_chile2010(tmp_path: Path, save: bool, variant: str):
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
    topo = _get_topography(runner)

    # Setup data for test
    runner.set_data()
    # runner.rundata.clawdata.num_output_times = 1
    runner.rundata.amrdata.amr_levels_max = 1

    if variant == "standard":
        runner.rundata.topo_data.topofiles = [[2, topo.path]]
    else:
        nc_path = tmp_path / f"chile_topo_{variant}.nc"
        meta = _write_netcdf_variant(topo, variant, nc_path)
        # Pass TopoMetadata as the third element so data.py writes the
        # descriptor block Fortran needs to interpret coordinate conventions.
        runner.rundata.topo_data.topofiles = [[4, nc_path, meta]]

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

    runner.check_gauge(save=save, gauge_id=32412)

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
