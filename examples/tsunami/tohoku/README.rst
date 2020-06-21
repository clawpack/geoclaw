
.. _geoclaw_examples_tsunami_tohoku:

Tsunami arising offshore of Japan, 5:46AM March 11, 2010 (UTC)
==============================================================

This example uses the dtopo file `fujii.txydz`, which is downloaded from the 
script `maketopo.py`.  This can be run via::

    $ make topo

which also downloads a topo file for the ocean bathymetry.
This bathymetry originally came from the NOAA National Geophysical Data
Center (NGDC), now NCEI (see `Sources of tsunami data
<http://www.clawpack.org/tsunamidata.html>`__).

**Creating kml files to view on Google Earth**

Users have the option of creating  KML files which can loaded into Google Earth.  To create 
these files, make sure that the appropriate figure number is set in `setplots.py` and that 
`plotdata.kml = True`.  

To view the results in Google Earth, use the command

    $ open _plots/Tohoku_2011.kmz 

Assuming you have Google Earth installed, this should open the file in the Google Earth browser. 

You can `unzip` the `kmz` file if you want to extract individual `kml`
files.

See `<http://www.clawpack.org/googleearth_plotting.html>`_ for
documentation.

Version
-------

- Updated for v5.7.0 on 20 June 2020

