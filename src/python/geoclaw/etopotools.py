
"""

Tools to download ETOPO topography/bathymetry data from NCEI (formerly NGDC).
See https://www.ncei.noaa.gov/products/etopo-global-relief-model

Two entry points are provided:

- :func:`fetch_etopo` -- the recommended path.  Reads an ETOPO netCDF DEM
  (ETOPO 2022 by default) into a :class:`~clawpack.geoclaw.topotools.Topography`
  via :func:`clawpack.geoclaw.topotools.fetch_remote_topo`.  This is the
  consistently-available, near-best-available-data source.

- :func:`etopo1_download` -- legacy.  Downloads a topo_type 3 (ASCII) file from
  the old NGDC WCS-proxy endpoint.  That endpoint is legacy and often flaky;
  prefer :func:`fetch_etopo` (or
  :func:`clawpack.geoclaw.topotools.fetch_remote_topo`) instead.
"""


def fetch_etopo(name='etopo22_30sec', crop_extent=None, coarsen=1, buffer=0,
                align=None, nc_params={}, verbose=False):
    r"""Fetch an ETOPO netCDF DEM as a `Topography`.

    Thin convenience wrapper over
    :func:`clawpack.geoclaw.topotools.fetch_remote_topo` for the ETOPO netCDF
    nicknames in ``topotools.remote_topo_urls`` (e.g. the default
    ``'etopo22_30sec'`` = ETOPO 2022 30 arcsecond, or ``'etopo1'``).

    :Input:

     - *name* (str) - nickname (key of ``topotools.remote_topo_urls``) or a URL
       to an ETOPO netCDF file.  Default ``'etopo22_30sec'``.
     - *crop_extent* ([x1, x2, y1, y2] or None) - requested crop in domain
       coordinates; ``None`` reads the whole file.
     - *coarsen* (int) - coarsening factor (1 = none).
     - *buffer* (int) - points to keep outside the crop on each side.
     - *align* (tuple) - alignment when coarsening; see ``Topography.crop``.
     - *nc_params* (dict) - options forwarded to the ``topo_type=4`` reader
       (e.g. ``z_var``, ``assume_units``); see ``Topography.read``.
     - *verbose* (bool) - if True, print the resolved source.

    :Output:

     - a :class:`~clawpack.geoclaw.topotools.Topography` object.
    """

    from clawpack.geoclaw import topotools

    return topotools.fetch_remote_topo(name, crop_extent=crop_extent,
                                       coarsen=coarsen, buffer=buffer,
                                       align=align, nc_params=nc_params,
                                       verbose=verbose)


def etopo1_download(xlimits, ylimits, dx=0.0166666666667, dy=None,
        output_dir='.', file_name=None, force=False, verbose=True,
        return_topo=None):

    """
    Download etopo1 topography from NCEI and save as a topo_type 3 file, then
    return it as a `Topography` object.

    .. note::
        This uses the old NGDC WCS-proxy endpoint, which is legacy and often
        flaky.  For a consistently-available, modern netCDF source prefer
        :func:`fetch_etopo` or
        :func:`clawpack.geoclaw.topotools.fetch_remote_topo`.

    :Inputs:

    - *xlimits*: tuple (x1, x2) limits in longitude
      Must either have -180 <= x1 < x2 <= 180
           or 180 <= x1 < x2 <= 360
           or -360 <= x1 < x2 <= -180
      To download topo for a region spanning longitude 180, you must
      download two separate files, one on each side.

    - *ylimits*: tuple (y1, y2) limits in latitude
      Must have -90 <= y1 < y2 <= 90.

    - *dx*: resolution in x, default is 1./60. degree = 1 arcminute.
    - *dy*: resolution in y, default is dy = dx.
    - *output_dir*: directory to store file, default is '.'
    - *file_name*: name of file, default is constructed from xlimits,ylimits
    - *force*: if True, download even if the file already exists.
    - *verbose*: if True, print info from clawpack.clawutil.data.get_remote_file
    - *return_topo*: deprecated and ignored; a `Topography` is always returned.

    :Output:

    - a :class:`~clawpack.geoclaw.topotools.Topography` object read from the
      downloaded topo_type 3 file.

    Note: New NGDC format gives cell-registered values, so shift the
    values `xllcorner` and `yllcorner` to the specified corner.

    **To do:** Check whether it is possible to specify grid-registered
    values as implied at http://www.ngdc.noaa.gov/mgg/global/global.html

    The `nodata_value` line expected by GeoClaw is now missing,
    so add this in too.
    """

    from clawpack.geoclaw import topotools
    from clawpack.clawutil.data import get_remote_file
    import os
    import warnings
    from numpy import round

    if return_topo is not None:
        warnings.warn(
            "etopo1_download's return_topo argument is deprecated and ignored; "
            "a Topography object is now always returned.",
            DeprecationWarning, stacklevel=2)

    format = '&format=aaigrid'   # topo_type 3

    if dy is None:
        dy = dx

    arcminute = 1/60.
    if abs(dx-arcminute)>1e-8 or abs(dy-arcminute)>1e-8:
        print('*** Warning: data may not be properly subsampled at')
        print('*** resolutions other than 1 arcminute, dx=dy=1/60.')

    x1,x2 = xlimits
    y1,y2 = ylimits

    if file_name is None:
        file_name = 'etopo1_%i_%i_%i_%i_%imin.asc' \
            % (int(round(x1)),int(round(x2)),int(round(y1)),int(round(y2)),\
              int(round(60*dx)))

    if (x1>=180) and (x1<x2) and (x2<=360):
        longitude_shift = -360.
    elif (x1>=-360) and (x1<x2) and (x2<=-180):
        longitude_shift = 360.
    else:
        longitude_shift = 0.
    x1 = x1 + longitude_shift
    x2 = x2 + longitude_shift

    if (x1<-180) or (x1>=x2) or (x2>180):
        raise ValueError("Require -180 <= x1 < x2 <= 180 or 180 <= x1 < x2 <=360")
    if (y1<-90) or (y1>=y2) or (y2>90):
        raise ValueError("Require -90 <= y1 < y2 <= 90")

    bbox = '&bbox=%1.4f,%1.4f,%1.4f,%1.4f' % (x1,y1,x2,y2)
    res = '&resx=%1.12f&resy=%1.12f' % (dx,dy)
    url = 'http://maps.ngdc.noaa.gov/mapviewer-support/wcs-proxy/wcs.groovy' \
            + '?request=getcoverage&version=1.0.0&service=wcs' \
            + '&coverage=etopo1&CRS=EPSG:4326' \
            + format + bbox + res

    file_path = os.path.join(output_dir,file_name)
    if os.path.exists(file_path) and (not force):
        print("Skipping download... file already exists: ",file_path)

    else:
        get_remote_file(url, output_dir=output_dir, file_name=file_name, \
                        verbose=verbose,force=force)

        x1 = x1 - longitude_shift   # shift back before writing header

        with open(file_path) as f:
            lines = f.readlines()
        if lines[2].split()[0] != 'xllcorner':
            print("*** Error downloading, check the file!")
        else:
            x1file = float(lines[2].split()[1])
            x2file = float(lines[3].split()[1])
            lines[2] = 'xllcorner    %1.12f\n' % x1file
            lines[3] = 'yllcorner    %1.12f\n' % x2file
            if 'nodata_value' not in lines[5]:
                lines = lines[:5] + ['nodata_value    -99999\n'] + lines[5:]
                print("Added nodata_value line")
            with open(file_path, 'w') as f:
                f.writelines(lines)
        print("Created file: ",file_path)

    topo = topotools.Topography()
    topo.read(file_path, topo_type=3)
    return topo
