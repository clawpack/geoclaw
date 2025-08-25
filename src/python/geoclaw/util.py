#!/usr/bin/env python

r"""
GeoClaw util Module  `$CLAW/geoclaw/src/python/geoclaw/util.py`

Module provides provides utility functions.


:Functions:

 - dms2decimal - Convert (degrees, minutes, seconds) to decimal degrees
 - dist_meters2latlong - Convert dx, dy distance in meters to degrees
 - dist_latlong2meters - Convert dx, dy distance in degrees to meters
 - haversine - Calculate the haversine based great circle distance
 - inv_haversine - Inverts the haversine distance
 - bearing - Compute the bearing from on location to another
 - gctransect - Compute a set of points on the great circle between two points
 - fetch_noaa_tide_data - Fetches water levels and tide predictions
"""

import io
import os
import os.path

import numpy
from urllib.parse import urlencode
from urllib.request import urlopen

# ==============================================================================
#  Constants
# ==============================================================================
from clawpack.geoclaw.data import Rearth, DEG2RAD, RAD2DEG, LAT2METER

NOAA_API_URL = 'https://tidesandcurrents.noaa.gov/api/datagetter'

# ==============================================================================
#  Functions for calculating Distances
# ==============================================================================
def dms2decimal(d,m,s,coord='N'):
    r"""Convert coordinates in (degrees, minutes, seconds) to decimal form.  
    
    If coord == 'S' or coord == 'W' then value is negated too.

    :Example: 

        >>> topotools.dms2decimal(7,30,36,'W')
        -7.51

    (Note that you might want to add 360 to resulting W coordinate
    if using E coordinates everywhere in a computation spanning date line.)

    :Returns: float

    """

    deg = d + m / 60.0 + s / 3600.0
    if coord in ['S','W']:
        deg = -deg

    return deg


def dist_latlong2meters(dx, dy, latitude=0.0):
    """Convert distance from degrees longitude-latitude to meters.

    Takes the distance described by *dx* and *dy* in degrees and converts it into
    distances in meters.

    returns (float, float) 

    """

    dym = Rearth * DEG2RAD * dy
    dxm = Rearth * numpy.cos(latitude * DEG2RAD) * dx * DEG2RAD

    return dxm,dym


def dist_meters2latlong(dx, dy, latitude=0.0):
    """Convert distance from meters to degrees of longitude-latitude.

    Takes the distance described by *dx* and *dy* in meters and converts it into
    distances in the longitudinal and latitudinal directions in degrees.  

    returns (float, float)

    """

    dxd = dx / (Rearth * numpy.cos(latitude * DEG2RAD)) * RAD2DEG
    dyd = dy * RAD2DEG / Rearth

    return dxd, dyd


def haversine(x0, y0, x1=None, y1=None, units='degrees'):

    """
    x0,y0 is assumed to be a point (or an array with the same shapes as x1,y1)
    x1,y1 is a point or two arrays of points (of the same dimension)
    returns array with same shape as x1 and y1 containing distance of each point
    from (x0,y0).

    For backward compatibility, also allows x0,y0 to be 2-tuples specifying
    two points, but this is not suggested since the notation is not consistent.
    """
    
    if x1 is None:
        # for backward compatibility, assume in this case that x0 and y0 
        # are tuples for the two desired points:
        assert len(x0)==len(y0)==2, "*** Unexpected input"
        x1,y1 = y0
        x0,y0 = x0

    if units == 'degrees':
        # convert to radians:
        x0 = x0*DEG2RAD
        y0 = y0*DEG2RAD
        x1 = x1*DEG2RAD
        y1 = y1*DEG2RAD

    dx = x1 - x0
    dy = y1 - y0

    # angle subtended by two points, using Haversine formula:
    dsigma = 2.0 * numpy.arcsin( numpy.sqrt( numpy.sin(0.5 * dy)**2   \
            + numpy.cos(y0) * numpy.cos(y1) * numpy.sin(0.5 * dx)**2))

    return Rearth * dsigma


def inv_haversine(d,x1,y1,y2,Rsphere=Rearth,units='degrees'):
    r"""Invert the Haversine function to find dx given a distance and point.


    Invert the haversine function to find dx given distance d and (x1,y1) and y2.
    The corresponding x2 can be x1+dx or x1-dx.
    May return NaN if no solution.
    """

    if units=='degrees':
        # convert to radians:
        x1 = x1 * RAD2DEG
        y1 = y1 * RAD2DEG
        y2 = y2 * RAD2DEG
    elif units != 'radians':
        raise Exception("unrecognized units")
    dsigma = d / Rsphere
    cos_dsigma = (numpy.cos(dsigma) - numpy.sin(y1)*numpy.sin(y2)) / (numpy.cos(y1)*numpy.cos(y2))
    dx = numpy.arccos(cos_dsigma)
    if units=='degrees':
        dx = dx * RAD2DEG
    return dx


def bearing(x0, y0, x1, y1, units='degrees', bearing_units='degrees'):

    """
    Compute the bearing from (x0,y0) to (x1,y1), i.e., the angle clockwise from
    due North of the great circle path from point 0 to 1.  

    The value returned is thus between 0 and 360 if bearing_units='degrees',
    or between 0 and 2*pi if bearing_units='radians'.

    Note: If using this to initialize a radially-symmetric 2d velocity on the
    sphere based on a radial velocity U(r), symmetric about (x0, y0), set:
        # lat-long assumed to be in degrees, r in meters
        r = haversine(x0,y0,x,y)
        beta = bearing(x0,y0,x,y,bearing_units='radians')
        u = U(r) * sin(beta)  # beta measured from North!
        v = U(r) * cos(beta)
    """
    from math import atan2, degrees, radians

    if units == 'degrees':
        # convert to radians:
        x0 = x0*DEG2RAD
        y0 = y0*DEG2RAD
        x1 = x1*DEG2RAD
        y1 = y1*DEG2RAD
    elif units != 'radians':
        raise Exception("unrecognized units")

    dx = x1 - x0
    dy = y1 - y0
    xx = numpy.cos(y0)*numpy.sin(y1) - numpy.sin(y0)*numpy.cos(y1)*numpy.cos(dx)
    yy = numpy.sin(dx) * numpy.cos(y1)
    b = atan2(yy, xx)   # in radians from North (between -pi and pi)

    beta = (degrees(b) + 360) % 360  # convert to degrees clockwise from North

    if bearing_units == 'radians':
        beta = radians(beta)  # convert to radians clockwise (0 to 2*pi)

    elif bearing_units != 'degrees':
        raise Exception("unrecognized bearing_units")

    return beta
    
def gctransect(x1,y1,x2,y2,npts,coords='W',units='degrees',Rearth=Rearth):
    r"""
    Given (longitude,latitude) pairs (x1,y1), (x2,y2) and npts
    Compute (x,y) for npts equally spaced points on the great circle connecting
    them.
    
    If coords='W' the points will all have -2*pi < x <= 0.
    If coords='E' the points will all have -0 <= x < 2*pi.
    With continuity at the date line x = \pm pi.

    Sample usage for 50 points on great circle from Tohoku to Crescent City:

        from clawpack.geoclaw import util,kmltools
        xtrans,ytrans = util.gctransect(142,37,-124.2,41.74,50,'W')
        kmltools.transect2kml(xtrans,ytrans)
        # then open transect.kml to view on Google Earth

    Based in part on formulas in
         https://math.stackexchange.com/questions/1783746

    """
    from numpy import pi,sin,cos,arcsin,sqrt,arctan2,array,dot,zeros,linspace
    
    assert coords in ['E','W'], '*** coords must be E or W'
    assert units in ['radians','degrees'], '*** units must be degrees or radians'
    
    d2r = pi/180
    if units=='degrees':
        # convert to radians:
        x1 = x1 * d2r
        y1 = y1 * d2r
        x2 = x2 * d2r
        y2 = y2 * d2r    

    # compute Cartesian coordinates in 3D:
    V1 = array([cos(x1)*cos(y1), sin(x1)*cos(y1), sin(y1)])
    V2 = array([cos(x2)*cos(y2), sin(x2)*cos(y2), sin(y2)])
    d = dot(V1,V2)
    beta = sqrt(1/(1-d**2))
    alpha = -beta*d
    W = alpha*V1 + beta*V2
    yW = arcsin(W[2])
    xW = arcsin(W[1]/cos(yW))

    # compute transect points:
    t2 = arcsin(1/beta)
    if d<0: t2 = pi - t2
    t = linspace(0, t2, npts)

    xtrans = zeros(npts)
    ytrans = zeros(npts)

    for j,tj in enumerate(t):
        V = cos(tj)*V1 + sin(tj)*W
        yV = arcsin(V[2])        # latitude between -pi/2 and pi/2
        xV = arctan2(V[1],V[0])  # longitude between -pi and pi

        if coords == 'W' and V[1]>0:
            # want -2*pi < xV <= 0, continuous at -pi
            xV = xV - 2*pi
        if coords == 'E' and V[1]<0:
            # want 0 <= xV < 2*pi, continuous at pi
            xV = xV + 2*pi    
                    
        if units=='degrees':
            xV = xV/d2r
            yV = yV/d2r
        xtrans[j] = xV
        ytrans[j] = yV
        
    return xtrans, ytrans


# ==============================================================================
#  Data Fetching and Mangling
# ==============================================================================
def get_netcdf_names(path, lookup_type='dim', user_mapping=None, verbose=False):
    r"""Determine names of dimensions/variables in the NetCDF file at *path*

    Note that if a file has multiple matching variables that the first it finds
    will be the one used.  It is best practice to provide the label if there
    are multiple available, e.g. 'msl' (mean sealevel pressure) and 'sp'
    (surface pressure) for pressure.

    Default dimension mapping:
     - "x": ["x", "longitude", "lon"]
     - "y": ["y", "latitude", "lat"]
     - "z": ["z"]
     - "t": ["t", "time"]

    Default variable mapping:
     - "wind_u":   ["wind_x", "wind_u", "u10", "u"]
     - "wind_v":   ["wind_y", "wind_v", "v10", "v"]
     - "pressure": ["pressure", "msl", "sp"]
     - "topo":     ["topo", "elevation", "z", "band1"]

    :Input:
     - path (Path)
     - lookup_type (str)
     - user_mapping (dict)
     - verbose (bool)

    :Output:
     - (dict) 
    """

    import xarray as xr

    # Open file and construct set of names
    data = xr.open_dataset(path)
    if 'dim' in lookup_type.lower():
        _var_mapping = {"x": (False, ["x", "longitude", "lon"]), 
                        "y": (False, ["y", "latitude", "lat"]),
                        "z": (False, ["z"]),
                        "t": (False, ["t", "time"])
                       }
        var_names = [var_name.lower() for var_name in data.dims]
        # Assume that we only want ('x'), weird
        if len(var_names) == 1:
            _var_mapping.pop(['y', 'z', 't'])
        # Assume that we only want ('x', 'y')
        elif len(var_names) == 2:
            _var_mapping.pop('z', 't')
        # Assume that we have ('x', 'y', 't')
        elif len(var_names) == 3:
            _var_mapping.pop('z')
        # Maybe a third or extraneous dimension is present?
        elif len(var_names) == 4:
            pass
        else:
            raise ValueError(f"{len(var_names)} present, not sure what to do.")
    elif 'var' in lookup_type.lower():
        _var_mapping = {"wind_u":   (False, ["wind_u", "wind_x", "u10", "u"]),
                        "wind_v":   (False, ["wind_v", "wind_y", "v10", "v"]),
                        "pressure": (False, ["pressure", "msl", "sp"]),
                        "topo":     (False, ["topo", "elevation", "z", "band1"])
               }
        var_names = [var_name.lower() for var_name in data.keys()]
    else:
        raise ValueError(f"Unknown lookup type {lookup_type}")

    # Use defaults here, checking to make sure they are in the dimension names
    if user_mapping:
        for (coord, names) in user_mapping.items():
            if isinstance(names, str):
                if names in var_names:
                    _var_mapping[coord] = (True, names)
            else:
                for (i, name) in enumerate(names):
                    if name in var_names:
                        _var_mapping[coord] = (True, name)
                        break

    # Search for names in default names
    for var_name in var_names:
        for (coord, values) in _var_mapping.items():
            if not values[0]:
                if var_name.lower() in values[1]:
                    _var_mapping[coord] = (True, var_name)
                    break

    # Check dims are matched, Cannot check for variable names
    if 'dim' in lookup_type.lower():
        for (coord, value) in _var_mapping.items():
            if not value[0]:
                raise ValueError(f"Could not find matching dim name in " +
                                 f"{path.name} with dimensions {var_names}.")

    # Construct final dictionary
    var_names = {}
    for (key, value) in _var_mapping.items():
        if value[0]:
            var_names[key] = value[-1]
    return var_names
    

def wrap_coords(input_path, output_path=None, force=False, dim_mapping=None):
    r"""Wrap the x coordinates of a NetCDF file from [0, 360] to [-180, 180]

    :Input:
     - *input_path* (path) Path to input file.
     - *output_path* (path) Path to output file.  Default is slightly modified
       version of *input_path* with 'wrap' appended.
     - *force* (bool) By default this routine will not overwrite a file that
       alreaady exists.  The value *force = True* will overwrite the file.  
       Default is False
     - *dim_mapping* (dict) Dimension mapping for all dimensions in the NetCDF
       file, not just 'x'.  Default is {}.
    """

    import xarray as xr

    if not output_path:
        output_path = (input_path.parent / 
                                f"{input_path.stem}_wrap{input_path.suffix}")

    if not output_path.exists() or force:
        dim_mapping = get_netcdf_names(input_path, lookup_type='dim', 
                                                       user_mapping=dim_mapping)
        ds = xr.open_dataset(input_path)
        lon_atrib = ds.coords[dim_mapping['x']].attrs
        ds.coords[dim_mapping['x']] = \
                                (ds.coords[dim_mapping['x']] + 180) % 360 - 180
        wrapped_ds = ds.sortby(ds[dim_mapping['x']])
        wrapped_ds.coords[dim_mapping['x']].attrs = lon_atrib
        wrapped_ds.to_netcdf(output_path)


def fetch_noaa_tide_data(station, begin_date, end_date, time_zone='GMT',
                         datum='STND', units='metric', cache_dir=None,
                         verbose=True):
    """Fetch water levels and tide predictions at given NOAA tide station.

    The data is returned in 6 minute intervals between the specified begin and
    end dates/times.  A complete specification of the NOAA CO-OPS API for Data
    Retrieval used to fetch the data can be found at:

        https://tidesandcurrents.noaa.gov/api/

    By default, retrieved data is cached in the geoclaw scratch directory
    located at:

        $CLAW/geoclaw/scratch

    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of retrieval
      - end_date (datetime): end of date/time range of retrieval

    :Optional Arguments:
      - time_zone (string): see NOAA API documentation for possible values
      - datum (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
      - cache_dir (string): alternative directory to use for caching data
      - verbose (bool): whether to output informational messages

    :Returns:
      - date_time (numpy.ndarray): times corresponding to retrieved data
      - water_level (numpy.ndarray): preliminary or verified water levels
      - prediction (numpy.ndarray): tide predictions
    """
    # use geoclaw scratch directory for caching by default
    if cache_dir is None:
        if 'CLAW' not in os.environ:
            raise ValueError('CLAW environment variable not set')
        claw_dir = os.environ['CLAW']
        cache_dir = os.path.join(claw_dir, 'geoclaw', 'scratch')

    def fetch(product, expected_header, col_idx, col_types):
        noaa_params = get_noaa_params(product)
        cache_path = get_cache_path(product)

        # use cached data if available
        if os.path.exists(cache_path):
            if verbose:
                print('Using cached {} data for station {}'.format(
                    product, station))
            return parse(cache_path, col_idx, col_types, header=True)

        # otherwise, retrieve data from NOAA and cache it
        if verbose:
            print('Fetching {} data from NOAA for station {}'.format(
                product, station))
        full_url = '{}?{}'.format(NOAA_API_URL, urlencode(noaa_params))
        with urlopen(full_url) as response:
            text = response.read().decode('utf-8')
            with io.StringIO(text) as data:
                # ensure that received header is correct
                header = data.readline().strip()
                if header != expected_header or 'Error' in text:
                    # if not, response contains error message
                    raise ValueError(text)

                # if there were no errors, then cache response
                save_to_cache(cache_path, text)

                return parse(data, col_idx, col_types, header=False)

    def get_noaa_params(product):
        noaa_date_fmt = '%Y%m%d %H:%M'
        noaa_params = {
            'product': product,
            'application': 'NOS.COOPS.TAC.WL',
            'format': 'csv',
            'station': station,
            'begin_date': begin_date.strftime(noaa_date_fmt),
            'end_date': end_date.strftime(noaa_date_fmt),
            'time_zone': time_zone,
            'datum': datum,
            'units': units
        }
        return noaa_params

    def get_cache_path(product):
        cache_date_fmt = '%Y%m%d%H%M'
        dates = '{}_{}'.format(begin_date.strftime(cache_date_fmt),
                               end_date.strftime(cache_date_fmt))
        filename = '{}_{}_{}'.format(time_zone, datum, units)
        abs_cache_dir = os.path.abspath(cache_dir)
        return os.path.join(abs_cache_dir, product, str(station), dates, 
                            filename)

    def save_to_cache(cache_path, data):
        # make parent directories if they do not exist
        parent_dir = os.path.dirname(cache_path)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

        # write data to cache file
        with open(cache_path, 'w') as cache_file:
            cache_file.write(data)

    def parse(data, col_idx, col_types, header):
        # read data into structured array, skipping header row if present
        a = numpy.genfromtxt(data, usecols=col_idx, dtype=col_types,
                             skip_header=int(header), delimiter=',',
                             missing_values='')

        # return tuple of columns
        return tuple(a[col] for col in a.dtype.names)

    # only need first two columns of data; first column contains date/time,
    # and second column contains corresponding value
    col_idx = (0, 1)
    col_types = 'datetime64[m], float'

    # fetch water levels and tide predictions
    try:
        date_time, water_level = fetch(
            'water_level', 'Date Time, Water Level, Sigma, O or I (for verified), F, R, L, Quality',
            col_idx, col_types)
    except:
        print('*** Fetching water_level failed, returning None')
        date_time = None
        water_level = None

    try:
        date_time2, prediction = fetch('predictions', 'Date Time, Prediction',
                                       col_idx, col_types)
        if date_time is None:
            date_time = date_time2

    except:
        print('*** Fetching prediction failed, returning None')
        date_time2 = None
        prediction = None

    # ensure that date/time ranges are the same
    if (date_time is not None) and (date_time2 is not None):
        if not numpy.array_equal(date_time, date_time2):
            raise ValueError('Received data for different times')

    return date_time, water_level, prediction
