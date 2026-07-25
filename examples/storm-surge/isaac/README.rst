
.. _geoclaw_examples_storm_surge_hurricane_isaac:

Storm Surge from Hurricane Isaac (parametric and gridded forcing)
=================================================================

This example contains the data and setup for running a storm surge forecast for
Hurricane Isaac (2012).  It doubles as a demonstration of GeoClaw's two
meteorological-forcing families driven from the *same* ATCF best-track file
(``bal092012.dat``).  The example can be run via::

    make all

which downloads the necessary topography, writes the storm forcing file, runs
the simulation, and plots the results.

Forcing variants
----------------

The forcing is selected in ``setrun.py`` through ``rundata.surge_data``.  By
default the example uses the **parametric** Holland 1980 model
(``storm_specification_type = 'holland80'``), which reads the ATCF track and
writes a compact GeoClaw storm file.

To use **gridded** forcing instead, set ``storm_specification_type = 'data'``
(equivalently, family ``"gridded"``).  Two gridded formats are supported and
shown in ``setrun.py``:

- **OWI / ASCII (NWS12)** -- a ``.PRE``/``.WIN`` pressure/wind file pair
  (``isaac.PRE`` / ``isaac.WIN``).  This is the default gridded branch.
- **NetCDF (NWS13)** -- a single CF-compliant ``.nc`` file of wind (``u10``,
  ``v10``) and mean-sea-level pressure (``msl``); see :ref:`netcdf_input` for
  the variable/coordinate conventions and automatic discovery.  Enable the
  commented NetCDF lines in ``setrun.py`` and build with NetCDF support
  (``-DNETCDF``).

The committed ``regression_data`` covers the ``holland80`` and OWI-ASCII
variants, and the test suite additionally exercises ERA5/NWS13 NetCDF forcing.

See :ref:`quick_surge` for a step-by-step surge walkthrough, :ref:`setrun_surge`
for the ``surge_data`` attribute reference, and :ref:`surgedata` for storm data
sources.
