try:
    import xarray as xr
except:
    'You must install xarray in order to use the xarray backends'
    raise

from clawpack.geoclaw.xarray_backends import FGOutBackend, FGMaxBackend

# epsg code for lat-lon
# Optionally, provide an epsg code to assign the associated coordinate system to the file.
# default behavior assigns no coordinate system. Note that if no coordinate system is provided,
# it will come in a GIS with coordinates associated with row and column number (not the x and y position
# encoded in the netcdf).

epsg_code = 4326

# An example of a fgout grid.
filename = '_output/fgout0001.b0001'
# provide the .bxxx file if binary format is used or the
# .qxxx file if ascii format is used.
# the format, fg number, and frame number are inferred from the filename.

ds = xr.open_dataset(filename, engine=FGOutBackend, backend_kwargs={'epsg':epsg_code})
# ds is now an xarray object. It can be interacted with directly or written to netcdf using
ds.to_netcdf('fgout0001_0001.nc')

# An example of a fgmax grid.
filename = "_output/fgmax0001.txt"
ds = xr.open_dataset(filename, engine=FGMaxBackend, backend_kwargs={'epsg':epsg_code})
ds.to_netcdf('fgmax0001.nc')
