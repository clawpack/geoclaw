#!/usr/bin/env python
# encoding: utf-8

r"""
GeoClaw topotools Module  `$CLAW/geoclaw/src/python/geoclaw/topotools.py`

Module provides several functions for reading, writing and manipulating
topography (bathymetry) files.

:Classes:
 - Topography

:Functions:

 - determine_topo_type
 - create_topo_func
 - topo1writer
 - topo2writer
 - topo3writer
 - swapheader


:TODO:
 - Add sub and super sampling capababilities
 - Add functions for creating topography based off a topo function, incorporate
   the create_topo_func into Topography class, maybe allow more broad
   initialization ability to the class to handle this?
 - Add more robust plotting capabilities
"""

import os

import numpy

import clawpack.geoclaw.util as util
import clawpack.clawutil.data
import clawpack.geoclaw.data

# ==============================================================================
#  Topography Related Functions
# ==============================================================================
# Recognized NetCDF attribute names that carry a vertical datum / reference.
_DATUM_ATTR_NAMES = ('vertical_datum', 'vertical_datum_name', 'datum',
                     'geospatial_vertical_datum')


def extract_datum(*attr_dicts):
    r"""Return the first recognized vertical-datum attribute, or None.

    Searches each mapping in *attr_dicts* (e.g. a NetCDF variable's attrs then
    the dataset's global attrs) for any of :data:`_DATUM_ATTR_NAMES` and
    returns the first value found.  Informational only; GeoClaw applies no
    vertical-datum transformation.
    """
    for attrs in attr_dicts:
        for name in _DATUM_ATTR_NAMES:
            if name in attrs:
                return attrs[name]
    return None


def determine_topo_type(path, default=None):
    r"""Using the file suffix of path, attempt to deterimine the topo type.

    :Input:

     - *path* (string) - Path to the file.  Can include archive extensions (they
       will be stripped off).
     - *default* (object) - Value returned if no suitable topo type was
       determined.  Default is *None*.

    returns integer between 1-3 or *default* if nothing matches.

    """

    extension = os.path.splitext(
                  clawpack.clawutil.data.strip_archive_extensions(path))[-1][1:]

    topo_type = default
    if extension[:2] == "tt" or extension[:8] == 'topotype':
        topo_type = int(extension[-1])
    elif extension == 'xyz':
        topo_type = 1
    elif extension == 'asc':
        topo_type = 3
    elif extension == 'txyz':
        topo_type = 1
    elif extension == 'nc':
        topo_type = 4

    return topo_type


def _netcdf_window_indices(coords, lo, hi, margin, n, stride=1):
    r"""Half-open index window ``[i0, i1)`` into 1-D monotonic *coords*.

    Used by :meth:`Topography.read` to push a ``crop_extent`` down to the
    NetCDF read so only the needed hyperslab is loaded from disk (rather than
    materializing a whole global variable and cropping afterward).

    :Input:
     - *coords* (ndarray) - 1-D coordinate array, monotonic ascending **or**
       descending (as stored in the file).
     - *lo*, *hi* (float) - requested inclusive coordinate bounds.
     - *margin* (int) - extra points to keep on each side so that a subsequent
       :meth:`Topography.crop` still has every point it needs (its own buffer
       plus the ``coarsen`` alignment search).
     - *n* (int) - length of *coords*.
     - *stride* (int) - read stride; the low index is snapped **down** to a
       multiple of *stride* so the strided sub-window shares the same phase as
       striding the full array, keeping the sampled grid identical.

    Returns ``(0, n)`` (the full range) when the interval does not overlap
    *coords*, mirroring ``crop()``'s fall-back of leaving the array uncropped
    when the filter region misses the topography.
    """
    # A monotonic array clipped to [lo, hi] yields a contiguous True block for
    # either sort order, so first/last True index bound the window.
    mask = (coords >= lo) & (coords <= hi)
    idx = numpy.nonzero(mask)[0]
    if idx.size == 0:
        return 0, n
    i0 = max(0, int(idx[0]) - margin)
    i1 = min(n, int(idx[-1]) + margin + 1)
    # Snap the low index down to a multiple of stride so coords[i0:i1:stride]
    # is a phase-aligned subset of coords[::stride].
    i0 -= i0 % stride
    if i1 <= i0:
        return 0, n
    return i0, i1


def create_topo_func(loc,verbose=False):
    """
    Given a 1-dimensional topography profile specfied by a set of (x,z)
    values, create a lambda function that when evaluated will give the
    topgraphy at the point (x,y).  (The resulting function is constant in y.)

    :Example:

        >>> f = create_topo_func(loc)
        >>> b = f(x, y)

    :Input:
     - *loc* (list) - Create a topography file with the profile denoted by the
       tuples inside of loc.  A sample set of points are shown below.  Note
       that the first value of the list is the x location and the second is
       the height of the topography. ::


        z (m)
        ^                                                  o loc[5]  o
        |
        |                                          loc[4]
        |--------------------------------------------o-----> x (m) (sea level)
        |
        |                                o loc[2] o loc[3]
        |
        |
        |                           o loc[1]
        |
        |
        |__________________o loc[0]
        0.0


    """

    cmd_str = "lambda x,y: (x <= %s) * %s" % (loc[0][0],loc[0][1])
    for i in range(0,len(loc)-1):
        loc_str = " + (%s < x) * (x <= %s)" % (loc[i][0],loc[i+1][0])
        loc_str = "".join((loc_str," * ((%s - %s) " % (loc[i][1],loc[i+1][1])))
        loc_str = "".join((loc_str," / (%s - %s)" % (loc[i][0],loc[i+1][0])))
        loc_str = "".join((loc_str," * (x - %s) + %s)" % (loc[i][0],loc[i][1])))
        cmd_str = "".join((cmd_str,loc_str))
    cmd_str = "".join((cmd_str," + (%s < x) * %s" % (loc[-1][0],loc[-1][1])))

    if verbose:
        print(cmd_str)
    return eval(cmd_str)


def topo1writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints):
    """
    Function topo1writer will write out the topofiles by evaluating the
    function topo on the grid specified by the other parameters.

    Assumes topo can be called on arrays X,Y produced by numpy.meshgrid.

    Output file is of "topotype1," which we use to refer to a file with
    (x,y,z) values on each line, progressing from upper left corner across
    rows, then down.
    """
    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)

    topography.write(outfile, topo_type=1)


def topo2writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
    r"""Write out a topo type 2 file by evaluating the function *topo*.

    This routine is here for backwards compatibility and simply creates a new
    topography object and writes it out.

    """

    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)

    topography.write(outfile, topo_type=2)


def topo3writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
    r"""Write out a topo type 3 file by evaluating the function *topo*.

    This routine is here for backwards compatibility and simply creates a new
    topography object and writes it out.

    """

    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)

    topography.write(outfile, topo_type=3)


def fetch_topo_url(url, local_fname=None, force=None, verbose=False,
                        ask_user=False):
    """
    DEPRECATED:  Use *clawpack.clawutil.data.get_remote_file* instead (see note below).

    Replaces get_topo function.

    Download a topo file from the web, provided the file does not
    already exist locally.

    :Input:
        - *url* (str) URL including file name
        - *local_fname* (str) name of local file to create.
          If *local_fname == None*, take file name from URL
        - *force* (bool) If False, prompt user before downloading.

    For GeoClaw examples, some topo files can be found in
    `http://www.geoclaw.org/topo`_
    See that website for a list of archived topo datasets.

    If force==False then prompt the user to make sure it's ok to download,

    If force==None then check for environment variable CLAW_TOPO_DOWNLOAD
    and if this exists use its value.  This is useful for the script
    python/run_examples.py that runs all examples so it won't stop to prompt.

    This routine has been deprecated in favor of
    *clawpack.clawutil.data.get_remote_file*.  All the functionality should be
    the same but calls the other routine internally.
    """

    if force is None:
        CTD = os.environ.get('CLAW_TOPO_DOWNLOAD', None)
        force = (CTD in [True, 'True'])

    if local_fname is not None:
        output_dir = os.path.dirname(local_fname)
        file_name = os.path.basename(local_fname)

    clawpack.clawutil.data.get_remote_file(url, output_dir=output_dir,
                                                file_name=file_name,
                                                force=force,
                                                verbose=verbose,
                                                ask_user=ask_user)


def get_topo(topo_fname, remote_directory, force=None):
    """
    DEPRECATED:  Use *clawpack.geoclaw.util.get_remote_file* instead

    Download a topo file from the web, provided the file does not
    already exist locally.

    remote_directory should be a URL.  For GeoClaw data it may be a
    subdirectory of  http://www.clawpack.org/geoclaw/topo
    See that website for a list of archived topo datasets.

    If force==False then prompt the user to make sure it's ok to download,
    with option to first get small file of metadata.

    If force==None then check for environment variable CLAW_TOPO_DOWNLOAD
    and if this exists use its value.  This is useful for the script
    python/run_examples.py that runs all examples so it won't stop to prompt.
    """

    url = remote_directory + '/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, force=force)


def swapheader(inputfile, outputfile):
    r"""Swap the order of key and value in header to value first.

    Note that this is a wrapper around functionality in the Topography class.

    """
    topo = Topography(inputfile)
    topo.write(outputfile)



# ==============================================================================
#  Topography class
# ==============================================================================
class Topography(object):

    r"""Base topography class.

    A class representing a single topography file.

    :Properties:

    Note: Modified to check the `grid_registration` when reading or writing
    topo files and properly deal with `llcorner` registration in which case
    the x,y data should be offset by dx/2, dy/2 from the lower left corner
    specified in the header of a DEM file.

    :Initialization:
         -

    :Examples:

        >>> import clawpack.geoclaw.topotools as topo
        >>> topo_file = topo.Topography()
        >>> topo_file.read('./topo.tt3', topo_type=3)
        >>> topo_file.plot()

    """

    @property
    def z(self):
        r"""A representation of the data as an 1d array."""
        if (self._z is None) and self.unstructured:
            self.read(mask=False)
        return self._z
    @z.setter
    def z(self, value):
        self._z = value
    @z.deleter
    def z(self):
        del self._z

    @property
    def Z(self):
        r"""A representation of the data as a 2d array."""
        if self._Z is None:
            self.generate_2d_topo(mask=False)
        return self._Z
    @Z.setter
    def Z(self, value):
        self._Z = value
    @Z.deleter
    def Z(self):
        del self._Z

    @property
    def x(self):
        r"""One dimensional coorindate array in x direction."""
        if self._x is None:
            self.read(mask=False)
        return self._x
    @x.setter
    def x(self, value):
        self._extent = None
        self._x = value
    @x.deleter
    def x(self):
        del self._x

    @property
    def X(self):
        r"""Two dimensional coordinate array in x direction."""
        if self._X is None:
            self.generate_2d_coordinates(mask=False)
        return self._X
    @X.setter
    def X(self, value):
        self._extent = None
        self._X = value
        self._x = numpy.nan
    @X.deleter
    def X(self):
        del self._X

    @property
    def y(self):
        r"""One dimensional coordinate array in y direction."""
        if self._y is None:
            self.read(mask=False)
        return self._y
    @y.setter
    def y(self, value):
        self._extent = None
        self._y = value
    @y.deleter
    def y(self):
        del self._y

    @property
    def Y(self):
        r"""Two dimensional coordinate array in y direction."""
        if self._Y is None:
            self.generate_2d_coordinates(mask=False)
        return self._Y
    @Y.setter
    def Y(self, value):
        self._extent = None
        self._Y = value
        self._y = numpy.nan
    @Y.deleter
    def Y(self):
        del self._Y

    @property
    def extent(self):
        r"""Extent of the topography."""
        if self._extent is None:
            self._extent = ( numpy.min(self.x), numpy.max(self.x),
                             numpy.min(self.y), numpy.max(self.y) )
        return self._extent
    @extent.setter
    def extent(self, value):
        self._extent = value

    @property
    def delta(self):
        r"""Spacing of data points."""
        if self._delta is None:
            if self.unstructured:

                # Calculate the smallest spacing between grid points
                dx = numpy.inf
                dy = numpy.inf
                num_comparisons = self.x.shape[0] - 1
                for i in range(self.x.shape[0]):
                    for j in range(num_comparisons):
                        dx = min(dx, numpy.abs(self.x[i + j + 1] - self.x[i]))
                        dy = min(dy, numpy.abs(self.y[i + j + 1] - self.y[i]))

                    num_comparisons -= 1
                self._delta = [dx, dy]
            else:
                # All other topography types should have equally spaced grid
                # points in each direction
                begin_delta = numpy.array([abs(self.x[1] - self.x[0]),
                                           abs(self.y[1] - self.y[0])])
                end_delta =   numpy.array([abs(self.x[-2] - self.x[-1]),
                                           abs(self.y[-2] - self.y[-1])])
                if not numpy.allclose(begin_delta, end_delta, 1e-8):
                    raise ValueError("Grid spacing delta not constant, ",
                                     "%s != %s." % (begin_delta, end_delta))

                dx = numpy.round(begin_delta[0], 15)
                dy = numpy.round(begin_delta[1], 15)
                self._delta = (dx, dy)
        return self._delta


    def __init__(self, path=None, topo_type=None, topo_func=None,
                       unstructured=False, **kwargs):
        r"""Topography initialization routine.

        See :class:`Topography` for more info.

        """

        super(Topography, self).__init__()

        self.path = path
        self.topo_func = topo_func
        self.topo_type = topo_type

        self.unstructured = unstructured
        # On-file / Fortran missing-data sentinel.  In memory, missing cells
        # are represented as NaN (see read()); this value is only used when
        # writing files and is read back from file headers.
        self.no_data_value = -99999

        # Optional vertical datum / reference metadata (e.g. 'MSL', 'NAVD88').
        # Informational only -- GeoClaw applies no vertical transformation; it
        # is populated from NetCDF attributes on read (when present) and
        # written back out on NetCDF write.
        self.datum = None

        # Data storage for only calculating array shapes when needed
        self._z = None
        self._Z = None
        self._x = None
        self._X = None
        self._y = None
        self._Y = None
        self._extent = None
        self._delta = None

        self.coordinate_transform = lambda x,y: (x,y)

        # Preprocessing attributes — applied by read() after data is loaded.
        #
        # BOUNDS TERMINOLOGY (consistent across Python and Fortran):
        #   * ``extent``      — the *actual* loaded-data bounds.  Read-only
        #                       @property (below), in DOMAIN coordinates.
        #   * ``crop_extent`` — the *requested* crop region, in DOMAIN
        #                       coordinates.  Mirrors Fortran ``tp_crop_extent``
        #                       and the storm module's ``met_crop_extent``.
        #   * ``crop_bounds`` — the same region in FILE coordinates, carried
        #                       only in the NetCDF (type-4) descriptor
        #                       (``netcdf_utils.FileMetadata.crop_bounds`` ->
        #                       Fortran ``nc_crop_bounds``).  The descriptor
        #                       writer converts crop_extent -> crop_bounds by
        #                       subtracting lon_wrap_offset/x_shift.
        # Convention: the ``_extent`` suffix is domain coords; ``_bounds`` is
        # file coords.  All 4-element vectors are ordered [x1, x2, y1, y2]
        # (x-pair then y-pair).  'crop_extent' is named to avoid shadowing the
        # 'extent' property above.
        # PATH NOTE: 'path' already exists as an instance attribute set above.
        # No topo_path alias is needed; callers should use self.path.
        self.crop_extent: list[float] | None = None  # [x1,x2,y1,y2]; None=full domain
        self.coarsen: int = 1
        self.buffer: float = 0.0
        self.align = None
        self.x_shift: float = 0.0
        self.y_shift: float = 0.0
        self.z_shift: float = 0.0
        self.negate_z: bool = False
        self._netcdf_meta = None  # TopoMetadata set by _normalize_topofiles for topo_entries() format

        if path:
            self.read(path=path, topo_type=topo_type,
                      unstructured=unstructured, **kwargs)

    def set_xyZ(self, X, Y, Z):
        r"""
        Set _x, _y, and _Z attributes and then generate X,Y,Z.

        If X,Y are 1d arrays, then shape of Z should be (len(Y), len(X)).

        Allow X,Y to be 2d arrays of shape Z.shape, in which case
        first extract x,y
        """

        if X.ndim == 1:
            x = X
        else:
            x = X[0,:]

        if Y.ndim == 1:
            y = Y
        else:
            y = Y[:,0]

        if Z.shape != (len(y),len(x)):
            raise ValueError("shape of Z should be (len(y), len(x))")

        diffx = numpy.diff(x)
        diffy = numpy.diff(y)
        dx = numpy.mean(diffx)
        dy = numpy.mean(diffy)
        if dy < 0:
            Y = numpy.flipud(Y)
            y = numpy.flipud(y)
            diffy = numpy.diff(y)
            dy = numpy.mean(diffy)
            Z = numpy.flipud(Z)
        if diffx.max()-diffx.min() > 1e-3*dx:
            print('diffx.max()-diffx.min() = ', diffx.max()-diffx.min())
            raise ValueError("x must be equally spaced for structured topo")
        if diffy.max()-diffy.min() > 1e-3*dy:
            print('diffy.max()-diffy.min() = ', diffy.max()-diffy.min())
            raise ValueError("y must be equally spaced for structured topo")

        self.unstructured = False
        self._x = x
        self._y = y
        self._Z = Z
        self._X = None
        self._Y = None
        self.generate_2d_coordinates()

        if X.ndim == 2:
            assert numpy.allclose(self.X, X), '*** X set incorrectly?'
        if Y.ndim == 2:
            assert numpy.allclose(self.Y, Y), '*** Y set incorrectly?'


    def generate_2d_topo(self, mask=False):
        r"""Generate a 2d array of the topo."""

        # Check to see if we need to generate these
        if self._Z is None:

            if self.unstructured:
                # Really no way to do this here with performing interpolation via
                # extract.  Note that if the interpolation is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d arrays, first interpolate the data and" \
                                 + " try to perform this operation again.")

            if self.path is not None:
                # RJL: why do we expect 1d z?
                if self._z is None:
                # Try to read the data, may not have done this yet
                    if self.topo_type is None:
                        self.topo_type = determine_topo_type(self.path)
                    if self.topo_type is None:
                        raise ValueError("topo_type must be specified")
                    self.read(path=self.path, topo_type=self.topo_type, mask=mask)
                    if self._Z is not None:
                        # We are done, the read function did our work
                        return

                # See if self._X and self._Y are already computed and use them if
                # available, otherwise just use self._x and self._y
                if self._X is not None and self._Y is not None:
                    new_shape = self._X.shape
                else:
                    new_shape = (self._x.shape[0], self._y.shape[0])
                # Reshape, note that the mask follows along with the new array
                self._Z = numpy.reshape(self._z, new_shape)

            elif self.topo_func is not None:
                # Generate topo via topo_func
                ## self._Z = numpy.flipud(self.topo_func(self.X, self.Y))
                ## RJL:  Don't flip -- leave so Z[i,j] has same dimensions as X,Y
                ## Othewise does not plot properly.
                self._Z = self.topo_func(self.X, self.Y)


    def generate_2d_coordinates(self, mask=False):
        r"""Generate 2d coordinate arrays."""

        # Check to see if we need to generate these
        if self._X is None and self._Y is None:

            # RJL: Added this to generate from _x and _y if available.
            # Correct?
            if (self._x is not None) and (self._y is not None):
                self._X,self._Y = numpy.meshgrid(self._x, self._y)

        if self._X is None and self._Y is None:
            if self.unstructured:
                # Really no way to do this here with performing interpolation via
                # extract.  Note that if the interpolation is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d coordinates, first interpolate the data" \
                                 + " and try to perform this operation again.")

            if self.path is not None:
                if abs(self.topo_type) == 1:
                    # Reading this topo_type should produce the X and Y arrays
                    self.read(mask=mask)
                elif abs(self.topo_type) in [2,3]:
                    if self._x is None or self._y is None:
                        # Try to read the data to get these, may not have been done yet
                        self.read(mask=mask)
                    # Generate arrays
                    self._X, self._Y = numpy.meshgrid(self._x, self._y)
                else:
                    raise ValueError("Unrecognized topo_type: %s" % self.topo_type)

            elif self.topo_func is not None:
                if self._x is None or self._y is None:
                    raise ValueError("The x and y arrays must be set to ",
                                     "create 2d coordinate arrays.")
                self._X, self._Y = numpy.meshgrid(self._x, self._y)


            # If masking has been requested try to get the mask first from
            # self._Z and then self._z
            if mask:
                if self._Z is None:
                    # Check to see if we really need to do anything here
                    if isinstance(self._z, numpy.ma.MaskedArray):
                        # Try to create self._Z
                        self.generate_2d_topo(mask=mask)

                if isinstance(self._Z, numpy.ma.MaskedArray):
                    # Use Z's mask for the X and Y coordinates
                    self._X = numpy.ma.MaskedArray(self._X, mask=self._Z.mask,
                                                                     copy=False)
                    self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Z.mask,
                                                                     copy=False)


    def read(self, path=None, topo_type=None, unstructured=False,
             mask=False, filter_region=None, force=False, stride=[1, 1],
             nc_params={}):
        r"""Read in the data from the object's *path* attribute.

        Stores the resulting data in one of the sets of *x*, *y*, and *z* or
        *X*, *Y*, and *Z*.

        :Input:
         - *path* (str)  file to read
         - *topo_type* (int) - GeoClaw format topo_type
         - *unstructured* (bool) - default is False for lat-long grids.
         - *mask* (bool) - whether to store as masked array for missing
           values (default if False)
         - *filter_region* (tuple)
         - *stride* (list) - List of strides for the x and y dimensions
           respectively.  Default is *[1, 1]*.  Note that this is only
           implemented for NetCDF reading currently.
         - *nc_params* (dict) - options for NetCDF (`topo_type=4`) reading:

             - `z_var` (str): name of the elevation variable, if it cannot be
               auto-detected by CF `standard_name` or common names.
             - `assume_units` (str): unit to assume for the elevation variable
               when the file has **no** `units` attribute (e.g. `"m"`).  Units
               are otherwise required and never silently assumed: a file whose
               elevation variable lacks `units`, or whose units are not meters,
               raises `ValueError` (GeoClaw does not convert on read; pre-
               convert non-meter data to meters first).

        The first three might have already been set when instatiating object.

        """

        if (path is None) and (self.path is None):
            raise ValueError("*** Need to set path for file to read")

        if path:
            self.path = path   # set or perhaps reset
            self.topo_type = None  # force resetting below

        if unstructured:
            self.unstructured = unstructured

        # Check if the path is a URL and fetch data if needed or forced
        #if "http" in self.path:
        #    fetch_topo_url(self.path)
        # RJL: should switch to util.get_remote_file, but after fetching
        # still need to read it in, which that routine does not do.
        # Do we really want to support this?  Seems better for user
        # to fetch and store as desired filename and then read file.


        if self.topo_type is None:
            if topo_type is not None:
                self.topo_type = topo_type
            else:
                # Try to look at suffix for type
                self.topo_type = determine_topo_type(self.path)
                if self.topo_type is None:
                    #self.topo_type = 3
                    raise ValueError("topo_type must be specified")

        if self.unstructured:
            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)
            points = []
            values = []

            # Filter region if requested
            if filter_region is not None:
                for coordinate in data:
                    if filter_region[0] <= coordinate[0] <= filter_region[1]:
                        if filter_region[2] <= coordinate[1] <= filter_region[3]:
                            points.append(coordinate[0:2])
                            values.append(coordinate[2])

                if len(points) == 0:
                    raise Exception("No points were found inside requested " \
                                  + "filter region.")

                # Cast lists as ndarrays
                self._x = numpy.array(points[:,0])
                self._y = numpy.array(points[:,1])
                self._z = numpy.array(values)

            else:
                self._x = data[:,0]
                self._y = data[:,1]
                self._z = data[:,2]

        else:
            # Data is in one of the GeoClaw supported formats
            if abs(self.topo_type) == 1:
                import warnings
                warnings.warn(
                    "topo_type=1 is deprecated. Convert to topo_type=2, 3, or 4:\n"
                    "  topo.read()  # load the type-1 file\n"
                    "  topo.write('output.tt2', topo_type=2)  # save as type 2\n"
                    "Note: topo_type=1 assumes regularly spaced data. Genuinely "
                    "unstructured (scattered) point data must be gridded externally "
                    "(e.g., scipy.interpolate, GMT) before use with GeoClaw.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                _preprocessing_requested = (
                    self.crop_extent is not None
                    or self.coarsen != 1
                    or self.buffer != 0.0
                    or self.align is not None
                    or self.x_shift != 0.0
                    or self.y_shift != 0.0
                    or self.z_shift != 0.0
                    or self.negate_z
                )
                if _preprocessing_requested:
                    raise NotImplementedError(
                        "Preprocessing attributes (crop_extent, coarsen, buffer, align, "
                        "shift) are not supported for topo_type=1. Convert to type 2/3/4 first:\n"
                        "  topo_raw = Topography(path=self.path, topo_type=1)\n"
                        "  topo_raw.read()\n"
                        "  topo_raw.write('converted.tt2', topo_type=2)"
                    )
                data = numpy.loadtxt(self.path)
                N = [0,0]
                y0 = data[0,1]
                for (n, y) in enumerate(data[1:,1]):
                    if y != y0:
                        N[1] = n + 1
                        break
                N[0] = data.shape[0] // N[1]

                self._x = data[:N[1],0]
                self._y = data[::N[1],1]
                self._Z = numpy.flipud(data[:,2].reshape(N))
                dx = self.X[0,1] - self.X[0,0]
                dy = self.Y[1,0] - self.Y[0,0]
                self._delta = (dx,dy)

            elif abs(self.topo_type) in [2,3]:
                # Get header information
                N = self.read_header()  # note this also sets self._extent
                                        # self._x, self._y, self._delta,
                                        # and  self.grid_registration

                if abs(self.topo_type) == 2:
                    # Data is read in as a single column, reshape it
                    self._Z = numpy.loadtxt(self.path, skiprows=6).reshape(N[1],N[0])
                    self._Z = numpy.flipud(self._Z)
                elif abs(self.topo_type) == 3:
                    # Data is read in starting at the top right corner
                    self._Z = numpy.flipud(numpy.loadtxt(self.path, skiprows=6))

            elif abs(self.topo_type) == 4:
                from clawpack.geoclaw import netcdf_utils as _ncutils
                from clawpack.geoclaw.units import (
                    convert as _units_convert,
                    GEOCLAW_NETCDF_UNITS as _NC_UNITS,
                )

                # Allow explicit variable name override via nc_params (backward
                # compat with old 'z_var' key).
                _var_hint: str | None = nc_params.get('z_var', None)
                # Opt-in escape hatch for a file with no 'units' attribute;
                # units are otherwise required, never silently assumed.
                _assume_units: str | None = nc_params.get('assume_units', None)
                # Opt-out of the post-conversion magnitude sanity check.
                _skip_sanity: bool = nc_params.get('skip_sanity_check', False)

                with _ncutils.TopoInspector(
                    self.path, var_name=_var_hint, assume_units=_assume_units,
                    skip_sanity_check=_skip_sanity,
                ) as inspector:
                    # Auto-detect elevation variable if not provided
                    if inspector.var_name is None:
                        inspector.var_name = inspector._find_topo_var_name()

                    # Get CF coordinate metadata (no fill-value data scan here;
                    # fill values are handled below via NaN replacement)
                    _meta = inspector.inspect(inspector.var_name)
                    # Check and record units; warns if conversion is needed.
                    # In-memory read converts recognized non-meter units to
                    # meters (the conversion block below applies the factor);
                    # the Fortran descriptor path applies the same factor via
                    # the descriptor scale_factor.
                    _source_units = inspector._check_topo_units()

                    ds = inspector.ds
                    _x_name = _meta.x_name
                    _y_name = _meta.y_name
                    _var_name = inspector.var_name

                    # Record an optional vertical datum from CF/common
                    # attributes (informational only; no transformation).
                    self.datum = extract_datum(ds[_var_name].attrs, ds.attrs)

                    # Load the 1-D coordinate arrays in full: these are
                    # O(nx)+O(ny) (a few MB even for a global grid), unlike the
                    # O(nx*ny) elevation variable.  Kept in file order for now;
                    # any N→S flip is applied below after windowing.
                    _lon_full = numpy.asarray(ds[_x_name].values, dtype=float)
                    _lat_full = numpy.asarray(ds[_y_name].values, dtype=float)

                    # Load variable, squeezing singleton non-spatial dims
                    _da = ds[_var_name]
                    for _dim in list(_da.dims):
                        if _dim not in (_x_name, _y_name):
                            if _da.sizes[_dim] == 1:
                                _da = _da.isel({_dim: 0})
                            else:
                                raise ValueError(
                                    f"NetCDF variable '{_var_name}' has "
                                    f"non-singleton dimension '{_dim}' "
                                    f"(size {_da.sizes[_dim]}).  Cannot load "
                                    f"as static topography."
                                )

                    # Transpose to (lat, lon) = (y, x) order expected by
                    # Topography.
                    _da = _da.transpose(_y_name, _x_name)

                    # Push both `stride` and `crop_extent` down to xarray's lazy
                    # indexing so the NetCDF backend reads ONLY the requested
                    # hyperslab.  Otherwise `_da.values` materializes the whole
                    # variable (a global DEM is many GB, and CF fill decoding
                    # promotes it to float), which is prohibitively slow and can
                    # exhaust memory even when the caller asked for a small
                    # subset.  The crop window is computed on the cheap 1-D
                    # coordinate arrays and expanded by a margin so the post-read
                    # crop() below still reproduces its exact bounds/buffer/
                    # align/coarsen result on the in-memory sub-window.
                    _nx = _lon_full.size
                    _ny = _lat_full.size
                    if self.crop_extent is not None:
                        _x1, _x2, _y1, _y2 = self.crop_extent
                        _margin = (int(self.buffer) + 1) * max(int(self.coarsen), 1)
                        _i0, _i1 = _netcdf_window_indices(
                            _lon_full, _x1, _x2, _margin, _nx, stride[0]
                        )
                        _j0, _j1 = _netcdf_window_indices(
                            _lat_full, _y1, _y2, _margin, _ny, stride[1]
                        )
                    else:
                        _i0, _i1 = 0, _nx
                        _j0, _j1 = 0, _ny

                    _da = _da.isel({
                        _y_name: slice(_j0, _j1, stride[1]),
                        _x_name: slice(_i0, _i1, stride[0]),
                    })
                    _z_vals = numpy.asarray(_da.values, dtype=float)
                    _lon_vals = _lon_full[_i0:_i1:stride[0]]
                    _lat_vals = _lat_full[_j0:_j1:stride[1]]

                    # Flip to S→N (y increasing) if file stores N→S
                    if not _meta.y_increasing:
                        _lat_vals = _lat_vals[::-1]
                        _z_vals = _z_vals[::-1, :]

                    # Apply unit conversion if source is not already meters
                    _contract = _NC_UNITS.get('topo', 'm')
                    _meters_aliases = frozenset(
                        {'m', 'meter', 'meters', 'metre', 'metres'}
                    )
                    if _source_units and _source_units not in _meters_aliases:
                        _canonical = _ncutils._normalize_cf_unit(_source_units)
                        if _canonical is not None:
                            _factor = _units_convert(1.0, _canonical, _contract)
                            _z_vals = _z_vals * _factor

                    # Magnitude sanity check on the resolved (meters) field.
                    if not _skip_sanity:
                        _ncutils._check_magnitude(
                            'topo',
                            float(numpy.nanmin(_z_vals)),
                            float(numpy.nanmax(_z_vals)),
                            var_name=_var_name, path=str(self.path),
                        )

                    # Decoded fill values are already NaN (from xarray
                    # mask_and_scale).  NaN is the in-memory missing-data
                    # representation, so they pass through unchanged.

                self._x = _lon_vals
                self._y = _lat_vals
                self._Z = _z_vals

            elif abs(self.topo_type) == 5:
                # GeoTIFF
                try:
                    import gdal
                except ImportError as e:
                    print("Reading GeoTIFF files requires GDAL.")
                    raise e

                data = gdal.Open(self.path)
                z = data.GetRasterBand(1).ReadAsArray()
                transform = data.GetGeoTransform()
                x_origin = transform[0]
                y_origin = transform[3]
                dx = transform[1]
                dy = -transform[5]

                self._Z = numpy.flipud(z)
                self._x = numpy.linspace(x_origin,
                                   x_origin + (z.shape[0] - 1) * dx, z.shape[0])
                self._y = numpy.linspace(y_origin - (z.shape[1] - 1) * dy,
                                   y_origin, z.shape[1])


            else:
                raise IOError("Unrecognized topo_type: %s" % self.topo_type)

            if self.topo_type < 0:
                # positive Z means distance below sea level for these
                # topo_type's, contrary to our convention, so negate:
                self._Z = -self._Z

            # Make sure these are set to None to force re-generating:
            self._X = None
            self._Y = None

            # Normalize missing data to NaN in memory.  The numeric
            # self.no_data_value is only the on-file/Fortran sentinel (written
            # back out by write() and read from file headers); in memory
            # missing cells are NaN so that min/max, arithmetic (e.g. z_shift)
            # and masking behave consistently across all topo_types.
            if not isinstance(self._Z, numpy.ma.MaskedArray):
                self._Z = numpy.where(self._Z == self.no_data_value,
                                      numpy.nan, self._Z)
            if mask:
                self._Z = numpy.ma.masked_invalid(self._Z)

            # Perform region filtering by delegating to crop() so the index
            # bounds are computed in exactly one place and are inclusive of the
            # filter_region edges (the previous inline slice dropped the upper
            # edge row/column).
            if filter_region is not None:
                _filtered = self.crop(filter_region=filter_region)
                if _filtered is not None:
                    self._x = _filtered._x
                    self._y = _filtered._y
                    self._Z = _filtered._Z
                    self._X = None
                    self._Y = None
                    self._extent = None
                    self._delta = None

            # ---------------------------------------------------------------
            # Apply preprocessing attributes in-memory (original file unchanged).
            # Fortran applies the same attributes independently in read_topo_file
            # and read_topo_settings so neither side needs to write a modified copy.
            # Order:
            #   1. negate_z
            #   2. z_shift   (missing cells are NaN and stay NaN under the shift)
            #   3. x_shift   (shift x array; Fortran shifts xlowtopo/xhitopo)
            #   3b. y_shift  (shift y array; Fortran shifts ylowtopo/yhitopo)
            #   4+5. crop + coarsen via self.crop() (Fortran: crop+buffer done,
            #        coarsen not yet implemented)
            # Steps are skipped when the attribute equals its default value.
            # ---------------------------------------------------------------
            if self.negate_z:
                self._Z = -self._Z
            if self.z_shift != 0.0:
                # Missing cells are NaN and remain NaN under the offset.
                self._Z = self._Z + self.z_shift
            if self.x_shift != 0.0:
                self._x = self._x + self.x_shift
                self._extent = None
            if self.y_shift != 0.0:
                self._y = self._y + self.y_shift
                self._extent = None
            if self.crop_extent is not None or self.coarsen > 1:
                _cropped = self.crop(
                    filter_region=self.crop_extent,
                    coarsen=int(self.coarsen),
                    buffer=int(self.buffer),
                    align=self.align,
                )
                if _cropped is not None:
                    self._x = _cropped._x
                    self._y = _cropped._y
                    self._Z = _cropped._Z
                    self._X = None
                    self._Y = None
                    self._extent = None
                    self._delta = None


    def read_header(self):
        r"""Read in header of topography file at path.

        If a value returns numpy.nan then the value was not retrievable.  Note
        that this routine can read in headers whose values and labels are
        swapped.

        """

        if abs(self.topo_type) in [2,3]:

            # Default values to track errors
            num_cells = [numpy.nan,numpy.nan]
            self._extent = [numpy.nan,numpy.nan,numpy.nan,numpy.nan]
            self._delta = numpy.nan

            with open(self.path, 'r') as topo_file:
                # Check to see if we need to flip the header values
                first_line = topo_file.readline()
                try:
                    num_cells[0] = int(first_line.split()[0])
                except ValueError:
                    # Assume the header is flipped from what we expect
                    num_cells[0] = int(first_line.split()[-1])
                    value_index = -1
                    label_index = 0
                else:
                    value_index = 0
                    label_index = -1

                num_cells[1] = int(topo_file.readline().split()[value_index])

                xline = topo_file.readline().split()
                xll = float(xline[value_index])
                # drop 'x' character and convert remaining string to lower case:
                x_registration = xline[label_index][1:].lower()

                yline = topo_file.readline().split()
                yll = float(yline[value_index])
                # drop 'y' character and convert remaining string to lower case:
                y_registration = yline[label_index][1:].lower()

                if x_registration == y_registration:
                    self.grid_registration = x_registration
                    # expect registration in ['llcorner', 'llcenter', 'lower']
                else:
                    raise IOError("x_registration and y_registration don't " \
                        + "match: %s,%s" % (x_registration, y_registration))

                # parse line allowing possibility of dx and dy (or just dx=dy)
                line = topo_file.readline()
                tokens = line.split()
                values = []
                for token in tokens:
                    try:
                        v = float(token)
                        values.append(v)
                    except:
                        pass
                dx = values[0]
                if len(values) == 1:
                    dy = dx   # only dx given
                elif len(values) == 2:
                    dy = values[1]
                    self._delta = (values[0], values[1])  # if dx,dy on line
                else:
                    raise IOError("Cannot parse dx,dy line: %s" % line)
                self._delta = (dx, dy)


                self.no_data_value = float(topo_file.readline().split()[value_index])

                x = numpy.linspace(xll, xll+(num_cells[0]-1)*dx, num_cells[0])
                y = numpy.linspace(yll, yll+(num_cells[1]-1)*dy, num_cells[1])
                if self.grid_registration in ['lower', 'llcenter']:
                    # x,y are cell center / data locations:
                    self._x = x
                    self._y = y
                elif self.grid_registration == 'llcorner':
                    # x,y are lower left corners:
                    # data points are offset by dx/2, dy/2
                    self._x = x + dx/2.
                    self._y = y + dy/2.
                    print('*** Note: since grid registration is llcorner,')
                    print('    will shift x,y values by (dx/2, dy/2) to cell centers')
                else:
                    # assume that x,y are cell center / data locations:
                    self._x = x
                    self._y = y
                    print('*** Warning: Unrecognized grid_registration: %s' \
                                    % self.grid_registration)
                    print('    Assuming x,y at grid points')

                # set extent based on data locations (not lower corner for 'llcorner')
                self._extent = [self._x[0],self._x[-1],self._y[0],self._y[-1]]

        elif abs(self.topo_type) == 4:
            # NetCDF: use NetCDFInspector for CF-aware coordinate detection.
            # Only coordinate arrays are loaded here; Z data is deferred.
            from clawpack.geoclaw import netcdf_utils as _ncutils

            with _ncutils.NetCDFInspector(self.path) as inspector:
                _x_name = inspector._find_x_name()
                _y_name = inspector._find_y_name()
                ds = inspector.ds

                _lon_vals = numpy.asarray(ds[_x_name].values, dtype=float)
                _lat_vals = numpy.asarray(ds[_y_name].values, dtype=float)

                # Normalise to S→N (increasing y) so extent/delta are consistent
                if not inspector._detect_y_increasing(_y_name):
                    _lat_vals = _lat_vals[::-1]

                # Record an optional vertical datum (informational) so it is
                # available without loading Z (e.g. for the consistency check
                # in TopographyData.write).
                self.datum = extract_datum(
                    *[ds[v].attrs for v in ds.data_vars], ds.attrs)

            self._x = _lon_vals
            self._y = _lat_vals
            self._extent = [self._x[0], self._x[-1], self._y[0], self._y[-1]]
            self._delta = (
                float(self._x[1] - self._x[0]),
                float(self._y[1] - self._y[0]),
            )
            num_cells = (len(self._x), len(self._y))

        elif abs(self.topo_type) == 5:
            # GeoTIFF
            try:
                import gdal
            except ImportError as e:
                print("Reading GeoTIFF files requires GDAL.")
                raise e

            data = gdal.Open(self.path)
            # z = data.GetRasterBand(1).ReadAsArray()
            transform = data.GetGeoTransform()
            x_origin = transform[0]
            y_origin = transform[3]
            dx = transform[1]
            dy = -transform[5]

            # self._Z = numpy.flipud(z)
            self._x = numpy.linspace(x_origin,
                               x_origin + (z.shape[0] - 1) * dx, z.shape[0])
            self._y = numpy.linspace(y_origin - (z.shape[0] - 1) * dy,
                               y_origin, z.shape[1])

        else:
            raise IOError("Cannot read header for topo_type %s" % self.topo_type)

        return num_cells

    def write(self, path, topo_type=None, no_data_value=None, fill_value=None,
                header_style='geoclaw', Z_format="%15.7e", grid_registration=None,
                z_dtype="float32", compression=None):
        r"""Write out a topography file to path of type *topo_type*.

        Writes out a topography file of topo type specified with *topo_type* or
        inferred from the output file's extension, defaulting to 3, to path
        from data in Z.  The rest of the arguments are used to write the header
        data.

        :Input:
         - *path* (str)  - file to write
         - *topo_type* (int) - GeoClaw format topo_type
           **Note:** this is second positional argument, agreeing with
           the read function in this class.  It was the third argument in
           GeoClaw version 5.3.1 and earlier.
         - *no_data_value* - values used to indicate missing data
         - *fill_value* (float) - value to use if filling a masked array
         - *header_style* (str) - indicates format of header lines
             'geoclaw' or 'default'  ==> write value then label
                                     with grid_registration == 'lower' as default
             'arcgis' or 'asc' ==> write label then value
                                   with grid_registration == 'llcorner' as default
                                   (needed for .asc files in ArcGIS)

         - *Z_format* (str) - string format to use for Z values
           The default format "%15.7e" gives at least millimeter precision
           for topography with abs(Z) < 10000 and results in
           smaller files than the previous default of "%22.15e" used in
           GeoClaw version 5.3.1 and earlier.  A shorter format can be used
           if the user knows there are fewer significant digits, e.g.
           etopo1 data is integers and so has a resolution of 1 meter.
           In this case a cropped or coarsened version might be written
           with `Z_format = "%7i"`, for example.
         - *grid_registration* (str) - 'lower', 'llcorner', 'llcenter'
                or None for defaults described above.
         - *z_dtype* (str) - on-disk dtype of the elevation variable when
                writing NetCDF (`topo_type=4`).  Default `"float32"`, which
                still gives sub-millimeter precision for Earth topography
                (abs(Z) < 10000 m) while halving file size; pass `"float64"`
                for full double precision.  Ignored for ASCII topo types.
                The elevation is always written with a CF `units = "m"`
                attribute (GeoClaw requires meters; see :ref:`topo_netcdf`).
         - *compression* - NetCDF (`topo_type=4`) zlib compression for the
                elevation variable.  `None`/`False` (default) writes
                uncompressed; `True` uses zlib level 1 + byte shuffle; an int
                selects the zlib complevel; a dict is passed through verbatim.
                The compressed file stays randomly readable and needs no reader
                change.  See `netcdf_utils.compression_encoding`.

        """

        # Determine topo type if not specified
        if topo_type is None:
            # Look at the suffix of the path and the object's topo_type
            # attribute to try to deterimine which to use, default to the path
            # version unless it did not work
            path_topo_type = determine_topo_type(path, default=-1)

            if self.topo_type is not None and path_topo_type == -1:
                topo_type = self.topo_type
            elif path_topo_type != -1:
                topo_type = path_topo_type
            else:
                # Default to 3 if all else fails
                topo_type = 3

        # Default to this object's no_data_value if the passed is None,
        # otherwise the argument will override the object's value or it will
        # default to -99999 (default for the class)
        if no_data_value is None:
            no_data_value = self.no_data_value

        # Check to see if masks have been applied to topography, if so
        # replace with fill_value (or  numpy.ma default value e.g. 1e+20)
        if isinstance(self.Z, numpy.ma.MaskedArray):
            if fill_value is not None:
                Z = self.Z.filled(fill_value)
            else:
                Z = self.Z.filled()
        else:
            Z = self.Z

        # check for NaNs:
        num_nan = numpy.isnan(Z).sum()
        if num_nan > 0:
            print('*** Z contains %i nan values, replacing with %s' \
                  % (num_nan, no_data_value))
            Z = numpy.where(numpy.isnan(Z), no_data_value, Z)

        # also fill self.z in the same way for unstructured?

        if self.unstructured:
            with open(path, 'w') as outfile:
                for (i, topo) in enumerate(self.z):
                    outfile.write("%s %s %s\n" % (self.x[i], self.y[i], topo))

        elif topo_type == 1:
            import warnings
            warnings.warn(
                "Writing topo_type=1 is deprecated. Prefer topo_type=2 or 3 for ASCII "
                "output, or topo_type=4 for NetCDF. Type-1 output will be removed in "
                "a future release.",
                DeprecationWarning,
                stacklevel=2,
            )
            with open(path, 'w') as outfile:
                for j in range(len(self.y)-1, -1, -1):
                    latitude = self.y[j]
                    for (i, longitude) in enumerate(self.x):
                        outfile.write("%s %s %s\n" % (longitude, latitude, self.Z[j,i]))

        elif topo_type == 2 or topo_type == 3:

            if grid_registration is None:
                if header_style in ['geoclaw','default']:
                    grid_registration = 'lower'
                elif header_style in ['arcgis','asc']:
                    grid_registration = 'llcorner'
                else:
                    raise ValueError("*** Unrecognized header_style")

            if grid_registration in ['lower','llcenter']:
                xlower = self.x[0]
                ylower = self.y[0]
            elif grid_registration == 'llcorner':
                xlower = self.x[0] - self.delta[0]/2.
                ylower = self.y[0] - self.delta[1]/2.
            else:
                raise ValueError('Unrecognized grid_registration: %s' \
                                % grid_registration)
            xlabel = 'x' + grid_registration
            ylabel = 'y' + grid_registration

            with open(path, 'w') as outfile:
                # Write out header
                if header_style in ['geoclaw','default']:
                    outfile.write('%6i                              ncols\n' % Z.shape[1])
                    outfile.write('%6i                              nrows\n' % Z.shape[0])
                    outfile.write('%22.15e              %s\n' % (xlower,xlabel))
                    outfile.write('%22.15e              %s\n' % (ylower,ylabel))
                    if abs(self.delta[0] - self.delta[1])/self.delta[0] < 1e-8:
                        # write only dx in usual case:
                        outfile.write('%22.15e              cellsize\n' \
                                % self.delta[0])
                    else:
                        # write both dx and dy if they differ:
                        outfile.write('%22.15e    %22.15e          cellsize\n' \
                                % (self.delta[0], self.delta[1]))
                    outfile.write('%10i                          nodata_value\n' % no_data_value)
                elif header_style in ['arcgis','asc']:
                    outfile.write('ncols  %6i\n' % Z.shape[1])
                    outfile.write('nrows  %6i\n' % Z.shape[0])
                    outfile.write('%s  %22.15e\n' % (xlabel,xlower))
                    outfile.write('%s  %22.15e\n' % (ylabel,ylower))
                    outfile.write('cellsize %22.15e\n'  % self.delta[0])
                    outfile.write('nodata_value  %10i\n' % no_data_value)
                else:
                    raise ValueError("*** Unrecognized header_style")

                # Write out topography data
                Z_flipped = numpy.flipud(Z)
                if topo_type == 2:
                    Z_format = Z_format + "\n"
                    for i in range(Z.shape[0]):
                        for j in range(Z.shape[1]):
                            outfile.write(Z_format % Z_flipped[i,j])
                elif topo_type == 3:
                    Z_format = Z_format + " "
                    for i in range(Z.shape[0]):
                        for j in range(Z.shape[1]):
                            outfile.write(Z_format % Z_flipped[i,j])
                        outfile.write("\n")
                del Z_flipped

        elif topo_type == 4:
            # Write a CF-compliant NetCDF file via xarray, normalized through
            # netcdf_utils.CFNormalizer so the output matches exactly what the
            # TopoInspector read path (and the Fortran descriptor) expect:
            # CF coordinate names (longitude/latitude) with standard_name/axis/
            # units, and a CF _FillValue rather than a custom no_data_value
            # attribute.  This is the writer counterpart to read(topo_type=4).
            import xarray as xr
            from clawpack.geoclaw.netcdf_utils import (CFNormalizer,
                                                       compression_encoding)

            # Coordinates are always stored as float64 (mirroring the dtopo
            # writer).  float32 lon/lat quantizes to decimeters-to-meters near
            # high magnitudes (e.g. ~1.7 m near 180 deg), which is a sizable
            # fraction of a 1/9" cell (~3.4 m) and exceeds a full cell at 1/27"
            # and finer -- adjacent points can collapse to the same value,
            # breaking the uniform-grid assumption.  float64 keeps coordinate
            # precision at the nanometer level regardless of resolution.
            elevation = xr.DataArray(
                Z,
                dims=("latitude", "longitude"),
                coords={"latitude": numpy.asarray(self.y, dtype=numpy.float64),
                        "longitude": numpy.asarray(self.x, dtype=numpy.float64)},
                name="elevation",
                attrs={
                    "standard_name": "height_above_reference_ellipsoid",
                    "long_name": "Elevation relative to sea level",
                    "units": "m",
                    "positive": "up",
                },
            )
            if self.datum is not None:
                elevation.attrs["vertical_datum"] = str(self.datum)
            ds = xr.Dataset({"elevation": elevation})
            ds.attrs.update({
                "Conventions": "CF-1.7",
                "title": "Topography Data",
                "institution": "Unknown",
                "source": "Unknown",
                "history": "",
                "references": "",
                "comment": "Created by GeoClaw",
            })

            # CFNormalizer adds standard_name/axis/units to the coordinate
            # variables and resolves fill-value attributes; running it here
            # guarantees the file is already in normalized form on read.
            ds = CFNormalizer(ds).normalize()

            # Encode the numeric file sentinel as the CF _FillValue so the
            # reader masks those cells (xarray decodes _FillValue -> NaN).
            # elevation is stored as *z_dtype* (float32 by default): topo values
            # are well under 10,000 m, so float32 still gives sub-millimeter
            # precision while halving file size; pass z_dtype='float64' for full
            # precision.  Coordinate variables must not carry xarray's default
            # NaN _FillValue (nonsensical for a monotonic axis), so it's
            # explicitly cleared for those too.
            coord_names = [name for name in ds.coords if name != "elevation"]
            elevation_encoding = {"_FillValue": no_data_value,
                                  "dtype": z_dtype}
            # Optional zlib compression; let netCDF auto-chunk the 2-D grid.
            elevation_encoding.update(compression_encoding(compression))
            encoding = {"elevation": elevation_encoding}
            encoding.update({name: {"_FillValue": None}
                             for name in coord_names})
            ds.to_netcdf(path, encoding=encoding)


        else:
            raise NotImplementedError("Output type %s not implemented." % topo_type)


    def plot(self, axes=None, contour_levels=None, contour_kwargs={},
             limits=None, cmap=None, add_colorbar=True,
             plot_box=False, long_lat=True, fig_kwargs={}, data_break=0.,
             cb_kwargs={}):
        r"""Plot the topography.

        :Input:
         - *axes* (matplotlib.pyplot.axes) - If passed in, plot will be
           added to this axes.  Otherwise a new plot figure will be created
           (using *fig_kwargs*) and a new *axes* object created and returned.
         - *contour_levels* (list) - levels for contour lines if these are
           to be added (default *None*).  Set to [0.] to plot shoreline.
         - *contour_kwargs* (dict) - keyword arguments to be passed to
           contour command, e.g. {'colors':'r', 'linestyles': '-'}.
           Default is empty dict.
         - *limits* (list) - (min, max) of topo values for color map.
           Defaults to None, in which case (self.Z.min(), self.Z.max()) used.
         - *cmap* (matplotlib.colors.Colormap) - colormap, defaults to
           specialized map with blues for bathymetry and green/browns for topo.
         - *fig_kwargs* (dict) - keyword arguments to be passed to figure.
         - *plot_box* (bool or color specifier) - If evaluates to True, plot
           a box around limits of this topo.
         - *long_lat* (bool) - If this is a longitude-latitude plot then set the
           aspect of the plot to compensate for stretching.  If not then the
           aspect is set to "equal".
         - *data_break* (float) - when default cmap is used, the value to use
           to break between water and land colormaps.
           Defaults to 0., but for some topo files may need to use e.g. 0.01
           Or may want to show plots at different tide stage.
         - *cb_kwargs* (dict) - keyword arguments to be passed to colorbar
           e.g. 'shrink', 'extend', 'label'.  Can also set 'title' for cbar

        :Output:
         - *axes* (matplotlib.pyplot.axes) - the axes on which plot created.

        Note that:
          - if *type(self.Z)* is *numpy.ma.MaskedArray* then *pcolor* is used,
          - if *type(self.Z)* is *numpy.ndarray* then *imshow* is used.
            (This is faster for large files)
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        import clawpack.visclaw.colormaps as colormaps
        from clawpack.visclaw import plottools

        # Create axes if needed
        if axes is None:
            fig = plt.figure(**fig_kwargs)
            axes = fig.add_subplot(111)

        # Turn off annoying offset
        axes.ticklabel_format(style="plain", useOffset=False)
        for label in axes.get_xticklabels():
            label.set_rotation(20)

        region_extent = self.extent


        if limits is None:
            if self.unstructured:
                topo_extent = (numpy.nanmin(self.z), numpy.nanmax(self.z))
            else:
                topo_extent = (numpy.nanmin(self.Z), numpy.nanmax(self.Z))
        else:
            topo_extent = limits

        # Create color map - assume shore is at z = data_break
        if cmap is None:
            land_cmap = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                                 0.25:[0.0,1.0,0.0],
                                                  0.5:[0.8,1.0,0.5],
                                                  1.0:[0.8,0.5,0.2]})
            sea_cmap = plt.get_cmap('Blues_r')
            if topo_extent[0] >= 0.0:
                cmap = land_cmap
                norm = colors.Normalize(vmin=0.0, vmax=topo_extent[1])
            elif topo_extent[1] <= 0.0:
                cmap = sea_cmap
                norm = colors.Normalize(vmin=topo_extent[0], vmax=0.0)
            else:
                cmap, norm = colormaps.add_colormaps((land_cmap, sea_cmap),
                                                     data_limits=topo_extent,
                                                     data_break=data_break)
        else:
            norm = colors.Normalize(vmin=topo_extent[0], vmax=topo_extent[1])


        if self.unstructured:
            plot = axes.scatter(self.x, self.y, c=self.z, cmap=cmap, norm=norm,
                                                marker=',', linewidths=(0.0,))
        else:
            plot = plottools.pcolorcells(self.X, self.Y, self.Z,
                                         ax=axes, norm=norm, cmap=cmap)
        if add_colorbar:
            try:
                # this kwarg can't be passed directly:
                cb_title = cb_kwargs.pop('title')
            except:
                cb_title = None

            cbar = plt.colorbar(plot, ax=axes, **cb_kwargs)

            if cb_title is not None:
                cbar.ax.set_title(cb_title)

            if 'label' not in cb_kwargs.keys():
                cbar.set_label('Topography (m)')

        # levels = range(0,int(-numpy.min(Z)),500)

        if (contour_levels is not None) and (not self.unstructured):
            axes.contour(self.X, self.Y, self.Z, levels=contour_levels,
                 **contour_kwargs)


        # expand extent to include full cells, which are centered at X,Y:
        x1 = self.x.min() - self.delta[0]/2.
        x2 = self.x.max() + self.delta[0]/2.
        y1 = self.y.min() - self.delta[1]/2.
        y2 = self.y.max() + self.delta[1]/2.

        axes.set_xlim(x1,x2)
        axes.set_ylim(y1,y2)

        if plot_box:
            # plot a box around this topography region
            if type(plot_box) is bool:
                color = 'm'
            else:
                # assume plot_box is a valid color:
                color = plot_box
            plt.plot([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], color=color)


        if long_lat:
            mean_lat = 0.5 * (region_extent[3] + region_extent[2])
            axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))
        else:
            axes.set_aspect('equal')

        return axes


    def interp_unstructured(self, fill_topo, extent=None, method='nearest',
                                   delta=None, delta_limit=20.0,
                                   no_data_value=-99999, buffer_length=100.0,
                                   proximity_radius=100.0,
                                   resolution_limit=2000):
        r"""Interpolate unstructured data on to regular grid.

        Function to interpolate the unstructured data in the topo object onto a
        structured grid.  Utilizes a bounding box plus a buffer of size
        *buffer_length* (meters) containing all data unless *extent is not None*
        is *True*.  Then uses the fill topography *fill_topo* to fill in the
        gaps in the unstructured data.  By default this is done by masking the
        fill data with the extents, the value *no_data_value* and if
        *proximity_radius* (meters) is not 0, by a radius of *proximity_radius*
        from all grid points in the object.  Stores the
        result in the *self.X*, *self.Y* and *self.Z* object attributes.  The
        resolution of the final grid is determined by calculating the minimum
        distance between all *x* and *y* data with a hard lower limit of
        *delta_limit* (meters).

        Note that the function *scipy.interpolate.griddata* does not respect
        masks so a call to *numpy.ma.MaskedArray.compressed()* must be made to
        remove the masked data.

        :Input:
         - *fill_topo* (list) - List of Topography objects to use as fill data
           in the projection.
         - *extent* (tuple) - A tuple defining the rectangle of the sub-section.
           Must be in the form (x lower,x upper,y lower, y upper).
         - *method* (string) - Method used for interpolation, valid methods are
           found in *scipy.interpolate.griddata*.  Default is *nearest*.
         - *delta* (tuple) - Directly set the grid spacing of the interpolation
           rather than determining it from the data itself.  Defaults to *None*
           which causes the method to determine this value itself.
           Should be a 2-tuple of floats (delta_x, delta_y).
         - *delta_limit* (float) - Limit of finest horizontal resolution,
           default is 20 meters.
         - *no_data_value* (float) - Value to use if no data was found to fill in a
           missing value, ignored if `method = 'nearest'`. Default is `-99999`.
         - *buffer_length* (float) - Buffer around bounding box, only applicable
           when *extent* is None.  Default is `100.0` meters.
         - *proximity_radius* (float) - Radius every unstructured data point
           used to mask the fill data with.  Default is `100.0` meters.
         - *resolution_limit* (int) - Limit the number of grid points in a
           single dimension.  Raises a *ValueError* if the limit is violated.
           Default value is `2000`.

        Sets this object's *unstructured* attribute to *False* if successful.

        """

        import scipy.interpolate as interpolate
        from scipy.spatial import cKDTree

        # Convert meter inputs to degrees
        mean_latitude = numpy.mean(self.y)
        buffer_degrees = util.dist_meters2latlong(buffer_length, 0.0,
                                                  mean_latitude)[0]
        delta_degrees = util.dist_meters2latlong(delta_limit, 0.0,
                                                 mean_latitude)[0]

        # Calculate new grid coordinates
        if extent is None:
            extent = [ numpy.min(self.x) - buffer_degrees,
                       numpy.max(self.x) + buffer_degrees,
                       numpy.min(self.y) - buffer_degrees,
                       numpy.max(self.y) + buffer_degrees ]
        if delta is None:
            delta_x = max( numpy.abs(self.x[1:] - self.x[:-1]).min(), delta_degrees)
            delta_y = max( numpy.abs(self.y[1:] - self.y[:-1]).min(), delta_degrees)
        else:
            try:
                delta_x, delta_y = delta   # tuple provided
            except TypeError:
                delta_x = delta_y = delta  # assume float provided

        N = ( numpy.ceil((extent[1] - extent[0]) / delta_x),
              numpy.ceil((extent[3] - extent[2]) / delta_y) )
        if not numpy.all(numpy.array(N) < resolution_limit):
            raise ValueError("Calculated resolution too high, N=%s!" % str(N))
        self._X, self._Y = numpy.meshgrid(
                                     numpy.linspace(extent[0], extent[1], int(N[0])),
                                     numpy.linspace(extent[2], extent[3], int(N[1])))

        # The object's own (real) unstructured points always survive.
        points = numpy.array([self.x, self.y]).transpose()
        values = numpy.asarray(self.z)

        # Build a KD-tree of the real points once for the proximity test: fill
        # points within proximity_radius of any real point are dropped so the
        # real data is preferred there.  The test is always against the
        # original points, not fill added by earlier fill_topo entries.
        if proximity_radius > 0.0:
            proximity_radius_deg = util.dist_meters2latlong(
                proximity_radius, 0.0, mean_latitude)[0]
            data_tree = cKDTree(points)

        # Mask each fill topography (structured or unstructured) and append the
        # surviving points.
        if not isinstance(fill_topo, list):
            fill_topo = [fill_topo]
        for topo in fill_topo:
            if topo.unstructured:
                x_fill, y_fill, z_fill = topo.x, topo.y, topo.z
            else:
                x_fill = topo.X.flatten()
                y_fill = topo.Y.flatten()
                z_fill = topo.Z.flatten()
            x_fill = numpy.asarray(x_fill)
            y_fill = numpy.asarray(y_fill)
            z_fill = numpy.asarray(z_fill)

            # Keep fill points inside the target extent with valid data
            # (missing data is NaN in memory; also honor a numeric
            # no_data_value).
            keep = ((x_fill >= extent[0]) & (x_fill <= extent[1]) &
                    (y_fill >= extent[2]) & (y_fill <= extent[3]) &
                    ~numpy.isnan(z_fill) & (z_fill != no_data_value))
            fill_points = numpy.column_stack((x_fill[keep], y_fill[keep]))
            fill_values = z_fill[keep]

            # Drop fill points within proximity_radius of a real data point.
            if proximity_radius > 0.0 and fill_points.shape[0] > 0:
                near = data_tree.query_ball_point(fill_points,
                                                  proximity_radius_deg)
                far = numpy.array([len(hits) == 0 for hits in near], dtype=bool)
                fill_points = fill_points[far]
                fill_values = fill_values[far]

            points = numpy.concatenate((fill_points, points))
            values = numpy.concatenate((fill_values, values))

        # Use specified interpolation
        self._Z = interpolate.griddata(points, values, (self.X, self.Y),
                                                                  method=method)

        self._extent = extent
        self._delta = (delta_x, delta_y)
        self.unstructured = False


    def in_poly(self, polygon):
        r"""Return a boolean mask of grid points inside *polygon*.

        :Input:
         - *polygon* - sequence of ``(x, y)`` vertices describing a closed
           polygon (the closing edge from the last vertex back to the first is
           implied).

        :Output:
         - *mask* (numpy.ndarray of bool) - ``True`` at points lying inside
           *polygon*, ``False`` elsewhere.  For a structured grid the mask has
           the same shape as ``self.X`` / ``self.Y``; for unstructured data it
           is 1-D over the scattered ``self.x`` / ``self.y`` points.

        Uses :class:`matplotlib.path.Path` for a robust point-in-polygon test:
        it handles concave polygons and is independent of vertex winding order.

        Example -- keep only topography inside a region of interest::

            mask = topo.in_poly(region_vertices)
            topo.Z[~mask] = numpy.nan
        """
        from matplotlib.path import Path

        path = Path(numpy.asarray(polygon, dtype=float))
        if self.unstructured:
            points = numpy.column_stack((self.x, self.y))
            return path.contains_points(points)
        points = numpy.column_stack((self.X.ravel(), self.Y.ravel()))
        return path.contains_points(points).reshape(self.X.shape)


    def replace_values(self, indices, value=numpy.nan, method='fill'):
        r"""Replace the Z values at *indices* in place using *method*.

        :Input:
         - *indices* - sequence of ``(i, j)`` index pairs identifying the
           cells to replace, e.g. the output of ``numpy.argwhere(condition)``.
         - *value* (float) - constant used when ``method == 'value'``.
           Default ``numpy.nan``.
         - *method* (str) - how to choose replacement values:

             - ``'value'``   - set the cells to the constant *value*.
             - ``'nearest'`` - nearest-neighbor value from the remaining
               (non-replaced) cells.
             - ``'linear'``  - linear interpolation from the remaining cells;
               cells outside the convex hull of the remaining data are left
               as ``numpy.nan``.
             - ``'fill'``    - replace each cell with the average of the
               nearest surrounding non-replaced cells, growing the search box
               until at least one is found (the default).

        ``'nearest'`` and ``'linear'`` interpolate in index space, which is
        equivalent to physical space for a regularly spaced grid.
        """
        indices = numpy.asarray(indices)
        if indices.size == 0:
            return
        bad = numpy.zeros(self.Z.shape, dtype=bool)
        bad[indices[:, 0], indices[:, 1]] = True

        if method == 'value':
            self.Z[bad] = value

        elif method in ('nearest', 'linear'):
            import scipy.interpolate as interpolate
            good = ~bad
            if not good.any():
                raise ValueError("Cannot interpolate: no valid data remains "
                                 "outside *indices*.")
            gi, gj = numpy.nonzero(good)
            bi, bj = numpy.nonzero(bad)
            self.Z[bad] = interpolate.griddata(
                numpy.column_stack((gi, gj)), self.Z[good],
                numpy.column_stack((bi, bj)), method=method)

        elif method == 'fill':
            # Average the nearest surrounding non-replaced cells, growing the
            # inf-norm search box until at least one good cell is found.
            ny, nx = self.Z.shape
            bad_pairs = set((int(i), int(j)) for i, j in indices)
            for i0, j0 in indices:
                i0, j0 = int(i0), int(j0)
                r = 0
                while r < max(ny, nx):
                    r += 1
                    summation = 0.0
                    num_points = 0
                    for i in range(max(0, i0 - r), min(i0 + r + 1, ny)):
                        for j in range(max(0, j0 - r), min(j0 + r + 1, nx)):
                            if (i, j) not in bad_pairs:
                                summation += self.Z[i, j]
                                num_points += 1
                    if num_points > 0:
                        self.Z[i0, j0] = summation / num_points
                        break

        else:
            raise ValueError("Unrecognized method %r; expected 'value', "
                             "'nearest', 'linear', or 'fill'." % method)


    def replace_no_data_values(self, value=numpy.nan, method='fill'):
        r"""Replace missing (NaN) cells in Z using *method*.

        Missing data is represented in memory as ``numpy.nan`` (the numeric
        ``no_data_value`` is only the on-file sentinel).  This locates those
        cells and replaces them via :meth:`replace_values`.

        :Input:
         - *value* (float) - constant used when ``method == 'value'``.
         - *method* (str) - one of ``'value'``, ``'nearest'``, ``'linear'``,
           or ``'fill'``; see :meth:`replace_values`.
        """
        no_data_indices = numpy.argwhere(numpy.isnan(self.Z))
        self.replace_values(no_data_indices, value=value, method=method)


    def smooth_data(self, indices, r=1):
        r"""Filter topo data at *indices* by averaging surrounding data.

        Surrounding data is considered within the ball of radius *r* in the
        inf-norm.  Acts as a low-band pass filter and removes oscillatory data.

        :Input:
         - *indices* (list)
         - *r* (int)

        :Output:
         None
        """

        index_range = [None, None]
        for index in indices:
            for n in range(2):
                index_range[n] = list(range(max(0, index[n] - r),
                                       min(index[n] + r + 1, self.Z.shape[n])))
            num_points = 0
            summation = 0.0
            for i in index_range[0]:
                for j in index_range[1]:
                    summation += self.Z[i,j]
                    num_points += 1
            if num_points > 0:
                self.Z[index[0], index[1]] = summation / num_points


    def crop(self, filter_region=None, coarsen=1, buffer=0, align=None):
        r"""Crop region to *filter_region*

        Create a new Topography object that is identical to this one but cropped
        to the region specified by filter_region

        :Input:

            - *filter_region* (tuple): (x1,x2,y1,y2) desired new extent
            - *coarsen* (int): coarsening factor (by subsampling)
            - *buffer* (int): when possible, have at least this many points
                                outside the filter_region on each side
            - *align* (tuple): (xalign,yalign) = desired alignment if coarsening

        Setting *buffer > 0* may be useful to insure that the
        computational domain lies entirely inside a cropped topo file
        (In GeoClaw, cell-centered topo values B are computed by integrating
        topo file values that are viewed as pointwise values, so topo
        point values are needed out to the domain edges.)

        When subsampling with *coarsen > 1*, the *align* parameter may be
        useful to insure that the subsampling starts at an appropriate
        index.  For example, if the original topo has

            topo.x = [0, 0.5, 1, 1.5, 2, 2.5]

        then coarsening by 2 would result in

            newtopo.x = [0, 1, 2]   # if align[0] is an integer

        or

            newtopo.x = [0.5, 1.5, 2.5]   # if align[0] is an integer + 0.5

        Often in GeoClaw, if the original topofile is aligned with
        integer longitudes and latitudes, for example, then we want
        any subsampled topo to have the same property.

        In general, it tries to choose a starting index so that

            (newtopo.x[0] - align[0]) / dx_new

        is an integer, where *dx_new* is the spacing of points in the
        new topo after coarsening.  This may not be possible, since it
        depends on the alignment of the original topography, in which
        case it will choose the index for which the misalignment is minimized.

        :TODO:
         - Currently this does not work for unstructured data, could in principle
         - This could be a special case of in_poly although that routine could
           leave the resulting topography as unstructured effectively.
        """

        if self.unstructured:
            raise NotImplementedError("*** Cannot currently crop unstructured topo")

        if filter_region is None:
            # only want to coarsen, so this is entire region:
            #filter_region = [self.x[0],self.x[-1],self.y[0],self.y[-1]]
            filter_region = self.extent

        xlower,xupper,ylower,yupper = filter_region

        dx,dy = self.delta
        dx_new = dx*coarsen
        dy_new = dy*coarsen

        # Find indices of topo arrays in filter_region:
        try:
            ilower = (self.x >= filter_region[0]).nonzero()[0][0]
            iupper = (self.x <= filter_region[1]).nonzero()[0][-1]
            jlower = (self.y >= filter_region[2]).nonzero()[0][0]
            jupper = (self.y <= filter_region[3]).nonzero()[0][-1]
        except:
            print('*** filter_region does not overlap topo')
            return None

        # shift indices if needed for alignment:
        if (coarsen > 1) and (align is not None):
            xs = numpy.array([self.x[ilower + i] for i in range(coarsen)])
            offsets = (xs - align[0]) / dx_new
            offsets_frac = offsets - numpy.round(offsets)
            ioffset = numpy.argmin(abs(offsets_frac))
            ilower = ilower + ioffset
            iupper = iupper - numpy.remainder(iupper-ilower, coarsen)
            #print(f'+++ shifted ilower by ioffset={ioffset} to {ilower}')

            ys = numpy.array([self.y[jlower + j] for j in range(coarsen)])
            offsets = (ys - align[1]) / dy_new
            offsets_frac = offsets - numpy.round(offsets)
            joffset = numpy.argmin(abs(offsets_frac))
            jlower = jlower + joffset
            jupper = jupper - numpy.remainder(jupper-jlower, coarsen)
            #print(f'+++ shifted jlower by joffset={joffset} to {jlower}')

        # buffer, checking limits of arrays:
        ilower = numpy.maximum(0, ilower - buffer*coarsen)
        jlower = numpy.maximum(0, jlower - buffer*coarsen)
        iupper = numpy.minimum(len(self.x)-1, iupper + buffer*coarsen) + 1
        jupper = numpy.minimum(len(self.y)-1, jupper + buffer*coarsen) + 1

        # Create new topography object:
        newtopo = Topography()

        newtopo._x = self._x[ilower:iupper:coarsen]
        newtopo._y = self._y[jlower:jupper:coarsen]

        # Force regeneration of 2d coordinate arrays and extent if needed
        newtopo._X = None
        newtopo._Y = None
        newtopo._extent = None

        # Modify Z array as well
        newtopo._Z = self._Z[jlower:jupper:coarsen, ilower:iupper:coarsen]

        newtopo.unstructured = self.unstructured
        newtopo.topo_type = self.topo_type

        # print "Cropped to %s by %s array"  % (len(newtopo.x),len(newtopo.y))

        if 0:
            # debugging checks:
            xlower_outside = (xlower - newtopo.x[0]) / dx_new
            ylower_outside = (ylower - newtopo.y[0]) / dy_new
            xupper_outside = (newtopo.x[-1] - xupper) / dx_new
            yupper_outside = (newtopo.y[-1] - yupper) / dy_new

            print(f'+++ fractions of cells outside should be between' \
                  + f' {buffer-1} and {buffer} since buffer={buffer}:')
            # note: the statement above is not true if filter_region extends
            # to or beyond the edges of the original topo self.extent
            print(f'+++ xlower_outside={xlower_outside},' \
                  + f' xupper_outside={xupper_outside}')
            print(f'+++ ylower_outside={ylower_outside},' \
                  + f' yupper_outside={yupper_outside}')

            if align is not None:
                xalign = (newtopo.x[0] - align[0])/dx_new
                yalign = (newtopo.y[0] - align[1])/dy_new
                print(f'+++ x alignment: {xalign} should be integer')
                print(f'+++ y alignment: {yalign} should be integer')

        return newtopo

    def make_shoreline_xy(self, sea_level=0):
        r"""
        Returns an array *shoreline_xy* with 2 columns containing x and y values
        for all segements of the shoreline (defined to be the contour
        where self.z = sea_level) separated by [nan,nan] pairs.
        This allows all shorelines to be quickly plotted via:

            >>> plot(shoreline_xy[:,0], shoreline_xy[:,1])

        The shoreline can be saved as a binary *.npy* file via:

            >>> numpy.save(filename, shoreline_xy)

        which is much smaller than the original topography file.
        Reload via:

            >>> shoreline_xy = numpy.load(filename)
        """

        import matplotlib.pyplot as plt

        x = self.x
        y = self.y
        Z = self.Z
        fig = plt.figure()
        c = plt.contour(x,y,Z,[sea_level])
        # c is the level 0 contour as list of arrays, one for each segement
        # catenate these together separated by array([nan,nan]):
        shoreline_xy = c.allsegs[0][0]  # first segment
        for k in range(1,len(c.allsegs[0])):
            shoreline_xy = numpy.vstack((shoreline_xy, \
                           numpy.array([numpy.nan,numpy.nan]), c.allsegs[0][k]))
        plt.close(fig)
        return shoreline_xy


    def make_function(self, interp_kwargs={}):
        """
        Create a function of (x,y) that returns the topo Z interpolated to a
        point (or to a 1D transect or 2D grid of points).

        :Inputs:
            *interp_kwargs*: dictionary of parameter values to be passed to
                             RegularGridInterpolator.  See defaults below.

        :Outputs:
            *topo_func*:  The function created

        See the docstring in topo_func below for details on what shapes
        its arguments (x,y) can be.

        """
        from scipy.interpolate import RegularGridInterpolator

        if 'method' not in interp_kwargs.keys():
            interp_kwargs['method'] = 'linear'
        if 'bounds_error' not in interp_kwargs.keys():
            interp_kwargs['bounds_error'] = False
        if 'fill_value' not in interp_kwargs.keys():
            interp_kwargs['fill_value'] = numpy.nan

        ZT = self.Z.T  # so indices refer to (x,y) rather than (y,x)
        topo_func1 = RegularGridInterpolator((self.x, self.y), ZT,
                                             **interp_kwargs)

        def topo_func(x,y):
            """
            Function that interpolates from topo to (x,y).  This function
            simplifies the function created by RegularGridInterpolator
            so the user can call topo_func(x,y) rather needing a tuple as
            input, and checks inputs for allowed shapes:

            x,y can both be scalars, in which case a scalar is returned.
            If one of x,y is a 1D array, the other can be:
                a 1D array of the same length or
                a scalar (which is equivalent to providing a 1D array of the
                          right length with the scalar value repeated).
            If both are 2D arrays, they should have the same shape.

            If at least one of x,y is an array, an array of the same shape
            is returned.
            """

            import numpy as np
            err_msg = f'*** unexpected combination of shapes for x and y'
            if np.isscalar(x):
                assert np.isscalar(y) or y.ndim == 1, err_msg
            elif np.isscalar(y):
                assert np.isscalar(x) or x.ndim == 1, err_msg
            else:
                assert x.shape == y.shape, err_msg
            return topo_func1((x,y))

        return topo_func



# Define convenience dictionary of URLs for some online DEMs in netCDF form:
remote_topo_urls = {}

# global 1 arcminute topography:
remote_topo_urls['etopo1'] = \
    'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO1_Ice_g_gmt4.nc'

# global 30 arcsecond topography from etopo 2022:
remote_topo_urls['etopo22_30sec'] = \
    'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_bed_elev_netcdf/ETOPO_2022_v1_30s_N90W180_bed.nc'

# some 1/3 arcsecond coastal modeling DEMs:
server = 'https://www.ngdc.noaa.gov/thredds/dodsC/regional/'
remote_topo_urls['astoria'] = server + 'astoria_13_mhw_2012.nc'
remote_topo_urls['puget_sound'] = server + 'puget_sound_13_mhw_2014.nc'
remote_topo_urls['port_townsend'] = server + 'port_townsend_13_mhw_2011.nc'
remote_topo_urls['strait_of_juan_de_fuca'] = \
    server + 'strait_of_juan_de_fuca_13_navd88_2015.nc'



def read_netcdf(path, zvar=None, extent='all', coarsen=1, return_topo=True,
                return_xarray=False, buffer=0, align=None, verbose=False):

    """
    :Input:

     - *path* (str) - Path to the file to read, or url to remote file,
       or a key into the topotools.remote_topo_urls dictionary.
     - *zvar* (str) - variable to read as Z=elevation.
       if None, will try 'Band1', 'z', 'elevation'.
     - *extent* - [x1,x2,y1,y2] for desired subset, or 'all' for entire file
     - *coarsen* (int) - factor to coarsen by, 1 by default.
     - *return_topo* (bool) - if True, return a topotools.Topography object.
       default is True
     - *return_xarray* (bool) - if True, return an xarray.Dataset object.
       default is False
     - *buffer* (int): when possible, have at least this many points
                        outside the filter_region on each side
     - *align* (tuple): (xalign,yalign) = desired alignment if coarsening
                        See the doc string for Topography.crop()

    :Output:
     - topo and/or xarray_ds depending on what was requested.
       (either a single object or a tuple of two objects.)

    If `return_xarray == True` then `xarray` is used to read the data,
    otherwise `netCDF4` is used directly.

    Sample usage:

        from clawpack.geoclaw import topotools
        extent = [-126,-122,46,49]
        path = 'etopo1'
        topo = topotools.read_netcdf(path, extent=extent, coarsen=2, \
                                     buffer=1, align=(-126,46), verbose=True)

        # results in topo.x = array([-126.03333333, -126., ...])

        # to plot:
        topo.plot()

        # to save topofile for input to GeoClaw:
        topo.write('etopo_sample_2min.tt3', topo_type=3, Z_format='%.0f')

    This should give a 2-minute resolution DEM of the Western Washington coast.
    Note that etopo1 Z values are integers (vertical resolution is 1 meter)
    and using `Z_format='%.0f'` will save as integers to minimize file size.

    Note that the newer etopo 2022 30 arcsecond DEM can be sampled using
    path = 'etopo22_30sec', but this topo is aligned differently with e.g.
    x = -126. falling half way between points. Also note that Z values in the
    newer dataset are no longer integers.
    """

    from numpy import array
    import netCDF4
    if return_xarray:
        import xarray

    # check if path is a key in the remote_topo_urls dictionary:
    if path in remote_topo_urls.keys():
        path = remote_topo_urls[path]

    if verbose:
        print("Will read netCDF data from \n    %s" % path)

    assert (type(coarsen) is int) and (coarsen >= 1), \
        '*** coarsen must be a positive integer'

    if return_xarray:
        f = xarray.open_dataset(path)
    else:
        f = netCDF4.Dataset(path, 'r')

    if 'lon' in f.variables:
        x = f.variables['lon']
    elif 'x' in f.variables:
        x = f.variables['x']
    else:
        print('*** f.variables = ',f.variables)
        raise ValueError("*** Unrecognized x, lon in netCDF file")

    if 'lat' in f.variables:
        y = f.variables['lat']
    elif 'y' in f.variables:
        y = f.variables['y']
    else:
        print('*** f.variables = ',f.variables)
        raise ValueError("*** Unrecognized y, lat in netCDF file")

    # for selecting subset based on extent, convert to arrays if netCDF4 used:
    #if not return_xarray:

    x = array(x)
    y = array(y)

    if zvar is None:
        if 'Band1' in f.variables:
            zvar = 'Band1'
        elif 'z' in f.variables:
            zvar = 'z'
        elif 'elevation' in f.variables:
            zvar = 'elevation'
        else:
            print('*** f.variables = ',f.variables)
            raise ValueError("*** Unrecognized zvar in netCDF file")


    if extent == 'all':
        ilower = 0
        iupper = len(x) - 1
        jlower = 0
        jupper = len(y) - 1
    else:
        x1,x2,y1,y2 = extent
        # find indices of x,y arrays for points lying within extent:
        iindex = numpy.where(numpy.logical_and(x >= x1, x <= x2))[0]
        jindex = numpy.where(numpy.logical_and(y >= y1, y <= y2))[0]
        ilower = iindex[0]
        iupper = iindex[-1]
        jlower = jindex[0]
        jupper = jindex[-1]

    dx_new = coarsen * (x[1] - x[0])
    dy_new = coarsen * (y[1] - y[0])

    # shift indices if needed for alignment:
    if (coarsen > 1) and (align is not None):
        xs = numpy.array([x[ilower + i] for i in range(coarsen)])
        offsets = (xs - align[0]) / dx_new
        offsets_frac = offsets - numpy.round(offsets)
        ioffset = numpy.argmin(abs(offsets_frac))
        ilower = ilower + ioffset
        iupper = iupper - numpy.remainder(iupper-ilower, coarsen)
        print(f'+++ shifted ilower by ioffset={ioffset} to {ilower}')

        ys = numpy.array([y[jlower + j] for j in range(coarsen)])
        offsets = (ys - align[1]) / dy_new
        offsets_frac = offsets - numpy.round(offsets)
        joffset = numpy.argmin(abs(offsets_frac))
        jlower = jlower + joffset
        jupper = jupper - numpy.remainder(jupper-jlower, coarsen)
        print(f'+++ shifted jlower by joffset={joffset} to {jlower}')

    # buffer, checking limits of arrays:
    i1 = numpy.maximum(0, ilower - buffer*coarsen)
    j1 = numpy.maximum(0, jlower - buffer*coarsen)
    i2 = numpy.minimum(len(x)-1, iupper + buffer*coarsen) + 1
    j2 = numpy.minimum(len(y)-1, jupper + buffer*coarsen) + 1

    xs = x[i1:i2:coarsen]
    ys = y[j1:j2:coarsen]
    Zs = f.variables[zvar][j1:j2:coarsen, i1:i2:coarsen]

    Zs = array(Zs)

    if 0:
        # debugging checks:
        xlower,xupper,ylower,yupper = extent
        dx_new = xs[1] - xs[0]
        dy_new = ys[1] - ys[0]
        xlower_outside = (xlower - xs[0]) / dx_new
        ylower_outside = (ylower - ys[0]) / dy_new
        xupper_outside = (xs[-1] - xupper) / dx_new
        yupper_outside = (ys[-1] - yupper) / dy_new

        print(f'+++ fractions of cells outside should be between' \
              + f' {buffer-1} and {buffer} since buffer={buffer}:')
        # note: the statement above is not true if filter_region extends
        # to or beyond the edges of the original topo self.extent
        print(f'+++ xlower_outside={xlower_outside},' \
              + f' xupper_outside={xupper_outside}')
        print(f'+++ ylower_outside={ylower_outside},' \
              + f' yupper_outside={yupper_outside}')

        if align is not None:
            xalign = (xs[0] - align[0])/dx_new
            yalign = (ys[0] - align[1])/dy_new
            print(f'+++ x alignment: {xalign} should be integer')
            print(f'+++ y alignment: {yalign} should be integer')

    if verbose:
        print('Returning a DEM with shape = %s' \
                % str(Zs.shape))
        print('x ranges from %.5f to %.5f with dx = %.8f' \
                % (xs[0], xs[-1], (xs[1]-xs[0])))
        print('y ranges from %.5f to %.5f with dy = %.8f' \
                % (ys[0], ys[-1], (ys[1]-ys[0])))
        if align is not None:
            xalign = (xs[0] - align[0])/dx_new
            yalign = (ys[0] - align[1])/dy_new
            print(f'aligned in x as requested if {xalign} is an integer')
            print(f'aligned in y as requested if {yalign} is an integer')
    output = None

    if return_topo:
        topo = Topography()
        topo.set_xyZ(xs,ys,Zs)
        output = topo

    if return_xarray:
        # Create a new xarray.Dataset with this subsampled, coarsened data:
        dims = (len(xs),len(ys))
        xarray_ds = xarray.Dataset({'z':(dims,Zs)}, coords={'lon':xs, 'lat':ys})
        if output is None:
            output = xarray_ds
        else:
            output = (topo, xarray_ds)

    return output
