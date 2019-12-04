
This test problem is based on examples/tsunami/bowl-slosh-netcdf.

This uses the `clawpack.geoclaw.topotools` capabilities to create the netCDF
topofile `bowl.nc`, using the Python module `netCDF4`.
This must be installed for this to work.
See `<https://unidata.github.io/netcdf4-python/netCDF4/index.html>`_.

Running the Fortran code also requires netCDF properly installed in
order to read in the `bowl.nc` file.  The `Makefile` refers to an
environment variable `NETCDF4_DIR` that must be set appropriately.
See `<https://www.unidata.ucar.edu/software/netcdf/docs-fortran/>`_.
