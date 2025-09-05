
.. _geoclaw_examples_tsunami_bowl-slosh-netcdf:

Sloshing water in a parabolic bowl with netCDF topofile
=======================================================

Waves in a parabolic bowl with a flat surface sloshing around.
An exact analytic solution is known in which the surface stays flat.

This is the same example as in examples/tsunami/bowl-slosh, but
using a topo file specified as a netCDF file rather than an ASCII raster.
See the README in that directory for more about the problem being solved.

To create the topo file before running the code::

    make topo

This uses the `clawpack.geoclaw.topotools` capabilities to create the netCDF
topofile `bowl.nc`, using the Python module `netCDF4`. 
This must be installed for this to work.
See `<https://unidata.github.io/netcdf4-python/netCDF4/index.html>`_.

Running the Fortran code also requires netCDF properly installed in 
order to read in the `bowl.nc` file.  The `Makefile` refers to an
environment variable `NETCDF4_DIR` that must be set appropriately.
See `<https://www.unidata.ucar.edu/software/netcdf/docs-fortran/>`_.


Version
-------

- Updated for v5.8.0
