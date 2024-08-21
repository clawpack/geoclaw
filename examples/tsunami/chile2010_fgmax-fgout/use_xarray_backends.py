import glob

try:
    import rioxarray
    import xarray as xr
except ImportError:
    "You must install xarray and rioxarray in order to use the xarray backends"
    raise

from clawpack.geoclaw.xarray_backends import FGMaxBackend, FGOutBackend

# epsg code for lat-lon
# Optionally, provide an epsg code to assign the associated coordinate system to the file.
# default behavior assigns no coordinate system. Note that if no coordinate system is provided,
# it will come in a GIS with coordinates associated with row and column number (not the x and y position
# encoded in the netcdf).

epsg_code = 4326

# An example of a fgout grid.
filename = "_output/fgout0001.b0001"
# provide the .bxxx file if binary format is used or the
# .qxxx file if ascii format is used.
# the format, fg number, and frame number are inferred from the filename.

ds = xr.open_dataset(
    filename,
    engine=FGOutBackend,
    backend_kwargs={
        "epsg": epsg_code,
        "qmap": "geoclaw",
        # qmap is the qmap specified to the fgout object in setrun.py see
        # the following documentation page for more details.
        # https://www.clawpack.org/dev/fgout.html#specifying-q-out-vars
        "dry_tolerance": None,
        # variables that are not eta and B are masked
        # where h>dry_tolerance. To turn this functionality
        # off set dry_tolerance = None.
    },
)
# ds is now an xarray object. It can be interacted with directly or written to netcdf using
ds.to_netcdf("fgout0001_0001.nc")

# It is possible to combine all fgout files into a single netcdf file
# using xr.open_mfdataset (requires dask) or xr.concat (does not require dask)
# https://docs.xarray.dev/en/stable/generated/xarray.open_mfdataset.html
# https://docs.xarray.dev/en/latest/generated/xarray.concat.html
# for instructions on installing xarray with dask see:
# https://docs.xarray.dev/en/latest/getting-started-guide/installing.html#instructions

fgout_files = glob.glob("_output/fgout0001.b*")

try:
    ds_all = xr.open_mfdataset(
        fgout_files,
        engine=FGOutBackend,
        backend_kwargs={"epsg": epsg_code, "qmap": "geoclaw"},
    )
except ValueError:  # if dask is not available, use xr.concat.
# if dask is not installed xr.open_mfdataset() will fail with something like
# ValueError: unrecognized chunk manager dask - must be one of: []
    fgouts = []
    for filename in fgout_files:
        ds = xr.open_dataset(
            filename,
            engine=FGOutBackend,
            backend_kwargs={"epsg": epsg_code, "qmap": "geoclaw"},
        )
        fgouts.append(ds)

    ds_all = xr.concat(fgouts, dim="time")

# save out.
ds_all.to_netcdf("fgout_all.nc")

# An example of a fgmax grid.
filename = "_output/fgmax0001.txt"
ds = xr.open_dataset(filename, engine=FGMaxBackend, backend_kwargs={"epsg": epsg_code})
ds.to_netcdf("fgmax0001.nc")

# To see the use of clipping, change the tfinal in setrun to something like 2*3600.0
# the fgmax0001_clipped.nc will only be the area where the wave arrived within the considered time.
filename = "_output/fgmax0001.txt"
ds = xr.open_dataset(
    filename, engine=FGMaxBackend, backend_kwargs={"epsg": epsg_code, "clip": True}
)
ds.to_netcdf("fgmax0001_clipped.nc")
