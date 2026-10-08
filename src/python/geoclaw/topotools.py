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
import warnings

import numpy

import clawpack.geoclaw.util as util
import clawpack.clawutil.data
import clawpack.geoclaw.data
from clawpack.geoclaw import coordinate_tools
from clawpack.geoclaw import gridded_input

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


def _crop_indices(x, y, crop_extent, coarsen, buffer, align):
    r"""Half-open index bounds for a crop+coarsen+align window.

    Returns ``(ilower, iupper, jlower, jupper)`` to be sliced as
    ``arr[jlower:jupper:coarsen, ilower:iupper:coarsen]`` (and the coordinate
    subsets ``x[ilower:iupper:coarsen]``, ``y[jlower:jupper:coarsen]``), or
    ``None`` if *crop_extent* does not overlap the data.

    This is the shared coarsen/align index arithmetic used by both
    :meth:`Topography.crop` (on the in-memory arrays) and the ``topo_type=4``
    read (on the cheap 1-D NetCDF coordinate arrays), so ASCII and NetCDF reads
    of the same data produce identical grids.

    :Input:
     - *x*, *y* (ndarray) - 1-D coordinate arrays, **ascending** (precondition;
       this routine contains no N->S/E->W flip logic).
     - *crop_extent* ([x1, x2, y1, y2]) - requested crop in the same coords.
     - *coarsen* (int) - subsampling factor; ``align`` only has effect when
       ``coarsen > 1``.
     - *buffer* (int) - grid points to keep outside *crop_extent* on each side
       (expanded by ``buffer*coarsen`` native points, as in :meth:`crop`).
     - *align* ((xalign, yalign) or None) - desired alignment when coarsening;
       ``None`` means no phase snap (start at the crop window).

    ``coarsen`` and ``buffer`` are assumed already ``int()``-coerced.

    :Raises:
     - *ValueError* if *crop_extent* is not increasing in both coordinates.  A
       descending longitude pair is the natural way to *spell* a crop across the
       antimeridian, and treating it as an ordinary window silently produced an
       empty grid, so it is rejected explicitly.
     - *ValueError* if *crop_extent* overlaps the data but contains no grid
       point (a window narrower than one cell, falling between two points).

    Warns when *crop_extent* extends beyond the data and is therefore clipped:
    a crop is never wrapped, only reduced.
    """
    # A descending pair is not an empty window; it is almost always an attempt
    # to cross the antimeridian.  Both index lookups below succeed when
    # crop_extent[0] > crop_extent[1], giving iupper < ilower and so a zero-size
    # slice -- an empty Topography whose `.extent` then raises something opaque
    # far from the cause.  Fail here instead.
    if crop_extent[0] >= crop_extent[1] or crop_extent[2] >= crop_extent[3]:
        raise ValueError(
            f"crop_extent must increase in both coordinates, got "
            f"{list(crop_extent)}. To cross the antimeridian use the "
            f"continuous spelling (e.g. [-211, -99]) rather than the wrapped "
            f"one ([170, -170]) -- and note that a Topography does not wrap on "
            f"its own: longitude wrapping is applied by the Fortran reader "
            f"from a descriptor written by TopoInspector.topo_entries(), which "
            f"emits two entries for a cross-seam crop. See the 'Region "
            f"terminology' note in the Topography docstring.")

    # dx/dy computed the same way as the `delta` property (round to 15 places),
    # so the align fractional-offset search matches crop() bit-for-bit.
    dx = numpy.round(abs(x[1] - x[0]), 15)
    dy = numpy.round(abs(y[1] - y[0]), 15)

    # Per-axis index math (crop window, align shift, buffer) lives in
    # coordinate_tools.crop_indices so the topo and dtopo crop paths share one
    # implementation; this function keeps only the two-axis vocabulary --
    # crop_extent validation, the clipping warning, and the (i, j) tuple.
    #
    # Both axes are resolved before either failure is reported, because the
    # original single try/except around all four index lookups let a
    # *non-overlap* on one axis win over a sub-cell window on the other (the
    # IndexError escaped first).  Raising per axis instead would turn that
    # None into a ValueError.
    _empty = None
    xwin = ywin = None
    for _axis, _coords, _lo, _hi, _delta, _align in (
            ('x', x, crop_extent[0], crop_extent[1], dx,
             None if align is None else align[0]),
            ('y', y, crop_extent[2], crop_extent[3], dy,
             None if align is None else align[1])):
        try:
            _win = coordinate_tools.crop_indices(
                _coords, _lo, _hi, _delta, coarsen, buffer, _align)
        except coordinate_tools.EmptyCropWindow as e:
            _empty = _empty or e
            _win = 'empty'
        if _axis == 'x':
            xwin = _win
        else:
            ywin = _win

    if (xwin is None) or (ywin is None):
        # crop_extent does not overlap the data.  Reported by the caller, which
        # knows whether the fall-back is "leave uncropped" or something else;
        # in particular no clipping warning here -- nothing was clipped.
        return None

    if _empty is not None:
        # The window overlaps the extent but falls strictly between two grid
        # points, so it contains no data.  Unlike a genuine non-overlap (which
        # has a documented full-file fallback) nothing pins this case, and
        # returning the full file for a crop the user narrowed *too far* is
        # never what was meant.  Re-raised here to name both axes and spacings.
        raise ValueError(
            f"crop_extent {list(crop_extent)} lies between grid points and "
            f"contains no data: the grid spacing is dx={dx}, dy={dy}. Widen "
            f"the crop to at least one cell, or use buffer= to include the "
            f"surrounding points.") from _empty

    # Silent clipping is the other half of the antimeridian confusion: a crop
    # written in continuous coordinates ([-211, -99] for a file on [-180, 180])
    # is not wrapped, it is quietly reduced to the part that exists.  Say so
    # when more than a cell is dropped; the result is still returned, because
    # over-wide crops are a legitimate and common way to say "all of this".
    _clipped = []
    if crop_extent[0] < x[0] - dx:
        _clipped.append(f"x1 {crop_extent[0]} -> {x[0]}")
    if crop_extent[1] > x[-1] + dx:
        _clipped.append(f"x2 {crop_extent[1]} -> {x[-1]}")
    if crop_extent[2] < y[0] - dy:
        _clipped.append(f"y1 {crop_extent[2]} -> {y[0]}")
    if crop_extent[3] > y[-1] + dy:
        _clipped.append(f"y2 {crop_extent[3]} -> {y[-1]}")
    if _clipped:
        warnings.warn(
            f"crop_extent {list(crop_extent)} extends past the data, which "
            f"covers x=[{x[0]}, {x[-1]}], y=[{y[0]}, {y[-1]}]; it was clipped "
            f"({', '.join(_clipped)}). A Topography does not wrap: to cross "
            f"the antimeridian use TopoInspector.topo_entries().")

    ilower, iupper = xwin
    jlower, jupper = ywin

    return int(ilower), int(iupper), int(jlower), int(jupper)


# The ascending-window -> file-order-slice mapping is format-neutral geometry
# shared by every NetCDF input reader, so it lives in coordinate_tools.
# Re-exported here under its original private name.
_axis_file_slice = coordinate_tools.axis_file_slice


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



# Sentinel for the deprecated crop-region kwargs (filter_region=, extent=).
# Their real default (None) is itself a valid user value, so a distinct sentinel
# is needed to tell "not passed" apart from "passed None".
_CROP_EXTENT_UNSET = object()

# Sentinel for read()'s align= kwarg.  align=None is a valid user value ("no
# phase snap"), and callers such as fetch_remote_topo set the self.align
# attribute *before* calling read(); a distinct sentinel lets read() tell "not
# passed" (leave self.align alone) apart from "passed None" (override it).
_ALIGN_UNSET = object()


def _resolve_crop_extent(crop_extent, deprecated):
    r"""Fold deprecated region kwargs onto ``crop_extent``.

    The requested-crop rectangle has one canonical name, ``crop_extent`` (see the
    "Region terminology" section of :class:`Topography`).  Older APIs spelled it
    ``filter_region`` or ``extent``; this helper maps those onto ``crop_extent``.

    :Input:
     - *crop_extent* - the value of the canonical ``crop_extent`` argument.
     - *deprecated* (dict) - maps each old kwarg name to the value it was called
       with, or ``_CROP_EXTENT_UNSET`` if it was not passed.

    For any old name that WAS passed, emit a ``DeprecationWarning`` and use its
    value as ``crop_extent``; raise ``TypeError`` if ``crop_extent`` is also
    supplied (ambiguous).
    """
    for name, value in deprecated.items():
        if value is _CROP_EXTENT_UNSET:
            continue
        if crop_extent is not None:
            raise TypeError(
                "Got both 'crop_extent' and the deprecated '%s'; "
                "pass only 'crop_extent'." % name)
        warnings.warn(
            "The '%s' argument is deprecated; use 'crop_extent' instead." % name,
            DeprecationWarning, stacklevel=3)
        crop_extent = value
    return crop_extent


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

    :Region terminology:

    Several attributes/arguments describe rectangular regions, all ordered
    ``[x1, x2, y1, y2]`` (x-pair then y-pair).  Two axes distinguish them:
    *role* (a derived result vs. a requested crop) and *coordinate frame*
    (domain vs. file).

    - ``extent`` -- the *actual* bounds of the loaded data, in DOMAIN
      coordinates.  A read-only, lazily-computed :func:`property` (a *result*,
      not an input); see :attr:`extent`.
    - ``crop_extent`` -- the *requested* crop rectangle, in DOMAIN coordinates.
      This is the single canonical name for the request: it is both a persisted
      attribute (default ``None`` = no crop) and the argument name accepted by
      :meth:`read`, :meth:`crop`, :meth:`interp_unstructured`, and
      :func:`fetch_remote_topo`.  ``read(crop_extent=r)`` is equivalent to
      setting the attribute and then reading.  Mirrors Fortran
      ``topo_crop_extent``.  (The older argument spellings ``filter_region`` and
      ``extent=`` are deprecated aliases.)
    - ``crop_bounds`` -- the same requested crop expressed in FILE coordinates.
      Used only in the NetCDF (type-4) layer
      (``netcdf_utils.FileMetadata.crop_bounds`` -> Fortran ``nc_crop_bounds``);
      converted from ``crop_extent`` by subtracting ``lon_wrap_offset``/
      ``x_shift``.
    - ``crop()`` -- the operation that turns a ``crop_extent`` request (plus
      ``coarsen``/``buffer``/``align``) into a new cropped object.

    Convention: the ``_extent`` suffix denotes domain coordinates; ``_bounds``
    denotes file coordinates.

    :The antimeridian, and what a cropped Topography represents:

    **A Topography never wraps.**  It holds one ascending ``x`` array and one
    ascending ``y`` array, so it can represent a rectangle in a continuous
    coordinate frame and nothing else.  Wrapping is not a property of the
    object; it is applied by the *Fortran* reader, from a ``lon_wrap_offset``
    that only exists in a NetCDF descriptor.  This is the distinction behind
    most antimeridian confusion, so it is worth being concrete about the three
    cases:

    1. **Ordinary ascending crop inside the file.**  The common case; nothing
       special happens.

    2. **Continuous spelling, e.g. ``crop_extent=[-211, -99, ...]`` for a file
       on ``[-180, 180]``.**  What this means depends on where it is used, and
       the two answers are different on purpose:

       - **Reading into a Topography** (``read``, ``crop``): the part that lies
         off the file is *clipped, not wrapped* -- you get ``[-180, -99]`` and
         a ``UserWarning`` saying so.  An in-memory ``Topography`` is one
         ascending array; it has nowhere to put the wrapped part.
       - **Writing to ``topo.data``** (``TopographyData.write``): the crop is
         *split across the seam*.  The writer emits one entry per side, each
         with its own ``lon_wrap_offset`` and file-coordinate ``crop_bounds``,
         and Fortran reassembles them into a single continuous region.

       So a cross-seam crop works from ordinary setrun code -- set
       ``crop_extent`` and append the ``Topography`` -- with no descriptor
       handling by the caller.  ``buffer``, ``coarsen``, ``align`` and the
       shifts are carried onto every entry.

    3. **Wrapped spelling, ``crop_extent=[170, -170, ...]``.**  Rejected
       wherever it appears, including at ``topo.data`` write time: Fortran
       would resolve ``crop_bounds = 170.0 -170.0`` to ``mx=0, my=0`` -- an
       empty topography, with no error.  Use the continuous spelling
       (``[-190, -170]``), which case 2 handles.

    :meth:`netcdf_utils.TopoInspector.topo_entries` is the underlying
    machinery, and is still available if you want the entries directly; the
    writer now calls it for you.  Latitude is never wrapped -- a crop whose
    latitude runs off the file is an error, since there is no seam to cross.

    The general shape of it: **Python is the single-rectangle case; the wrap
    lives in the Fortran interface.**  The same split explains why
    ``coordinate_system`` gates wrapping (a projected x axis in meters has no
    seam to cross) and why ``crop_bounds`` is in file coordinates while
    ``crop_extent`` is in domain coordinates.

    :Order of preprocessing operations:

    ``crop_extent`` -> ``align`` -> ``buffer`` -> ``coarsen``, applied in that
    order by both :meth:`crop` and the Fortran reader
    (``apply_align_buffer_coarsen``).  Because the strided subsample is last,
    **``buffer`` counts coarsened output points, not native file points**: the
    index window is widened by ``buffer * coarsen`` native points on each side,
    so ``buffer=2, coarsen=4`` adds 2 points to each edge of the result, not 8.

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
        r"""Actual bounds of the loaded data, ordered ``[x1, x2, y1, y2]``.

        This is a derived *result* (the min/max of the loaded ``x``/``y``), in
        domain coordinates -- not a requested crop; for the crop request see
        ``crop_extent`` and the "Region terminology" section of the class
        docstring.  Computed lazily and cached in ``_extent``; the cache is
        invalidated (set to ``None``) whenever ``x``/``y`` change or ``read()``
        reloads/crops the data.
        """
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
        # See the "Region terminology" section of the class docstring for the
        # extent / crop_extent / crop_bounds glossary and the Python<->Fortran
        # name mapping. Convention recap: the ``_extent`` suffix is domain
        # coords, ``_bounds`` is file coords, all ordered [x1, x2, y1, y2].
        # PATH NOTE: 'path' already exists as an instance attribute set above;
        # no topo_path alias is needed, callers should use self.path.
        self.crop_extent: list[float] | None = None  # [x1,x2,y1,y2]; None=full domain
        self.coarsen: int = 1
        self.buffer: int = 0
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
             mask=False, crop_extent=None, force=False,
             coarsen=None, align=_ALIGN_UNSET, buffer=None, stride=None,
             nc_params=None, filter_region=_CROP_EXTENT_UNSET):
        r"""Read in the data from the object's *path* attribute.

        Stores the resulting data in one of the sets of *x*, *y*, and *z* or
        *X*, *Y*, and *Z*.

        :Input:
         - *path* (str)  file to read
         - *topo_type* (int) - GeoClaw format topo_type
         - *unstructured* (bool) - default is False for lat-long grids.
         - *mask* (bool) - whether to store as masked array for missing
           values (default if False)
         - *crop_extent* ([x1, x2, y1, y2] or None) - requested crop region in
           domain coordinates (see the "Region terminology" section of the
           class docstring). Passing it here is equivalent to setting the
           ``crop_extent`` attribute before calling ``read()``; the crop is
           applied (together with ``coarsen``/``buffer``/``align``) via
           :meth:`crop`. Default ``None`` = no crop. The older ``filter_region``
           keyword is a deprecated alias.
         - *coarsen* (int) - subsampling factor (1 = no coarsening).  Applied
           identically for ASCII and NetCDF reads.  Passing it here is
           equivalent to setting the ``coarsen`` attribute before ``read()``.
           See :meth:`crop`.
         - *align* ((xalign, yalign) or None) - desired alignment when
           coarsening; see :meth:`crop`.  ``None`` (the default) means **no
           phase snap** -- subsampling starts at the crop window, matching
           ASCII/:meth:`crop`.  (This differs from the old NetCDF ``stride``
           behavior, which snapped to the file's grid origin.)  Pass e.g.
           ``align=[integer_lon, integer_lat]`` to lock the coarsened grid to a
           fixed lattice regardless of the requested ``crop_extent``.
         - *buffer* (int) - grid points to keep outside ``crop_extent`` on each
           side.  These are *output* points: with ``coarsen > 1`` the window
           grows by ``buffer * coarsen`` native points, so the result gains
           ``buffer`` points per edge, not ``buffer * coarsen``.  See
           :meth:`crop`.
         - *stride* (list or int) - **Deprecated**: use ``coarsen`` instead.
           A NetCDF-only knob that silently did nothing for ASCII reads and used
           a different alignment convention.  A scalar (or equal-valued list) is
           mapped onto ``coarsen``; per-axis striding is no longer supported.
         - *nc_params* (dict) - options for NetCDF (`topo_type=4`) reading:

             - `z_var` (str): name of the elevation variable, if it cannot be
               auto-detected by CF `standard_name` or common names.
             - `assume_units` (str): unit to assume for the elevation variable
               when the file has **no** `units` attribute (e.g. `"m"`), treated
               as if the file had declared it -- so `assume_units="km"` also
               converts.  Units are otherwise required and never silently
               assumed: a file whose elevation variable lacks `units` raises
               `ValueError`.  A *recognized* non-meter unit (e.g. `km`) is
               converted to meters on read, with a warning; an unrecognized
               unit raises.  See `dev/design/units_policy.md`.

        The first three might have already been set when instatiating object.

        """

        # None is the natural "no options" value and used to reach .get() as a
        # NoneType; it is also safer than a shared mutable default.
        if nc_params is None:
            nc_params = {}

        # A crop_extent passed here is equivalent to setting the attribute first;
        # fold the deprecated filter_region alias onto it, then store it so the
        # single attribute-driven crop below (and the type-4 pushdown) apply it.
        crop_extent = _resolve_crop_extent(crop_extent,
                                           {'filter_region': filter_region})
        if crop_extent is not None:
            self.crop_extent = crop_extent

        # Fold coarsen/align/buffer args onto the attributes (mirrors crop_extent
        # above): passing them to read() is equivalent to setting the attribute
        # first.  Sentinels distinguish "not passed" from an explicit value so a
        # caller that presets self.align/self.coarsen/self.buffer before read()
        # (e.g. fetch_remote_topo) is not silently clobbered.
        if coarsen is not None:
            self.coarsen = int(coarsen)
        if buffer is not None:
            self.buffer = int(buffer)
        if align is not _ALIGN_UNSET:
            self.align = align

        # `stride` is deprecated: a NetCDF-only knob that silently did nothing
        # for ASCII reads and used a different alignment convention than
        # crop()/coarsen.  Map it onto the unified scalar `coarsen`.
        if stride is not None:
            warnings.warn(
                "The 'stride' argument to Topography.read() is deprecated; use "
                "'coarsen' (a scalar subsampling factor) instead.  'coarsen' is "
                "applied identically for ASCII and NetCDF reads.",
                DeprecationWarning,
                stacklevel=2,
            )
            if numpy.ndim(stride) == 0:
                _stride = int(stride)
            else:
                _s = list(stride)
                if len(_s) == 0 or any(int(v) != int(_s[0]) for v in _s):
                    raise ValueError(
                        "Per-axis stride is no longer supported; 'coarsen' is a "
                        "single scalar factor.  Got stride=%r." % (stride,))
                _stride = int(_s[0])
            if _stride != 1:
                if self.coarsen != 1 and self.coarsen != _stride:
                    raise ValueError(
                        "Pass either 'stride' or 'coarsen', not both "
                        "(stride=%r, coarsen=%r)." % (stride, self.coarsen))
                self.coarsen = _stride

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
            # crop() already refuses unstructured data (NotImplementedError);
            # the read path used to attempt its own filter and then index a
            # Python list as if it were an array, dying with an unrelated
            # TypeError far from the cause.  Refuse consistently and early.
            _unstructured_preprocessing = (
                self.crop_extent is not None
                or self.coarsen != 1
                or self.buffer != 0
                or self.align is not None
            )
            if _unstructured_preprocessing:
                raise NotImplementedError(
                    "Preprocessing attributes (crop_extent, coarsen, buffer, "
                    "align) are not supported for unstructured data; "
                    "Topography.crop() refuses them for the same reason. Grid "
                    "the data first (e.g. interp_unstructured), then crop the "
                    "result.")

            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)
            points = []
            values = []

            # Filter region if requested (crop_extent, domain coords)
            if self.crop_extent is not None:
                for coordinate in data:
                    if self.crop_extent[0] <= coordinate[0] <= self.crop_extent[1]:
                        if self.crop_extent[2] <= coordinate[1] <= self.crop_extent[3]:
                            points.append(coordinate[0:2])
                            values.append(coordinate[2])

                if len(points) == 0:
                    raise Exception("No points were found inside requested " \
                                  + "crop_extent region.")

                # Cast lists as ndarrays
                self._x = numpy.array(points[:,0])
                self._y = numpy.array(points[:,1])
                self._z = numpy.array(values)

            else:
                self._x = data[:,0]
                self._y = data[:,1]
                self._z = data[:,2]

        else:
            # Data is in one of the GeoClaw supported formats.  Format-specific
            # reading is delegated to a registered GriddedReader (gridded_input);
            # the shared post-read processing below is identical for every
            # format.  Adding a format = registering a new adapter there.
            reader = gridded_input.get_reader(self.topo_type)
            _result = reader.read_window(self, nc_params=nc_params)
            self._x = _result.x
            self._y = _result.y
            self._Z = _result.Z
            if _result.delta is not None:
                self._delta = _result.delta
            if _result.datum is not gridded_input._UNSET:
                self.datum = _result.datum

            if self.topo_type < 0:
                # positive Z means distance below sea level for these
                # topo_type's, contrary to our convention, so negate:
                self._Z = -self._Z

            # Make sure these are set to None to force re-generating:
            self._X = None
            self._Y = None
            # _extent and _delta are derived from _x/_y, which were just
            # replaced.  read_header() populates them from the file header, so
            # without this a cropped topo_type=4 read (whose crop is applied
            # while reading the hyperslab, bypassing the property setters)
            # reported the *full file* extent alongside cropped data.
            self._extent = None
            self._delta = None

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
            #        coarsen not yet implemented).  SKIPPED for a reader that
            #        declares applies_preprocessing (the NetCDF reader applies
            #        crop+coarsen+align+buffer via _crop_indices while reading
            #        the hyperslab, so running crop() again would double-coarsen).
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
            if not reader.applies_preprocessing \
                    and (self.crop_extent is not None or self.coarsen > 1):
                _cropped = self.crop(
                    crop_extent=self.crop_extent,
                    coarsen=int(self.coarsen),
                    buffer=int(self.buffer),
                    align=self.align,
                )
                if _cropped is None:
                    # crop() already warned about the non-overlap; say what the
                    # consequence is here, because keeping the *full* file is
                    # surprising and the Fortran reader would instead abort.
                    warnings.warn(
                        f"crop_extent {list(self.crop_extent)} did not overlap "
                        f"{self.path}; the full file was read uncropped. The "
                        f"Fortran reader treats this as fatal.")
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


    def interp_unstructured(self, fill_topo, crop_extent=None, method='nearest',
                                   delta=None, delta_limit=20.0,
                                   no_data_value=-99999, buffer_length=100.0,
                                   proximity_radius=100.0,
                                   resolution_limit=2000,
                                   extent=_CROP_EXTENT_UNSET):
        r"""Interpolate unstructured data on to regular grid.

        Function to interpolate the unstructured data in the topo object onto a
        structured grid.  Utilizes a bounding box plus a buffer of size
        *buffer_length* (meters) containing all data unless *crop_extent is not
        None* is *True*.  Then uses the fill topography *fill_topo* to fill in the
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
         - *crop_extent* (tuple) - A tuple defining the rectangle of the
           sub-section, in the form (x1, x2, y1, y2). Default ``None`` uses the
           data bounding box plus *buffer_length*. The older ``extent`` keyword
           is a deprecated alias.
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
           when *crop_extent* is None.  Default is `100.0` meters.
         - *proximity_radius* (float) - Radius every unstructured data point
           used to mask the fill data with.  Default is `100.0` meters.
         - *resolution_limit* (int) - Limit the number of grid points in a
           single dimension.  Raises a *ValueError* if the limit is violated.
           Default value is `2000`.

        Sets this object's *unstructured* attribute to *False* if successful.

        """

        crop_extent = _resolve_crop_extent(crop_extent, {'extent': extent})
        # Internal working name for the interpolation output bounding box.
        extent = crop_extent

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


    def crop(self, crop_extent=None, coarsen=1, buffer=0, align=None,
             filter_region=_CROP_EXTENT_UNSET):
        r"""Crop region to *crop_extent*

        Create a new Topography object that is identical to this one but cropped
        to the region specified by *crop_extent* (see the "Region terminology"
        section of the class docstring).

        :Input:

            - *crop_extent* (tuple): (x1,x2,y1,y2) desired new extent, in domain
              coordinates. Default ``None`` crops to the current ``extent`` (so
              only *coarsen* has effect). The older ``filter_region`` keyword is
              a deprecated alias.
            - *coarsen* (int): coarsening factor (by subsampling). Truncated to
              an integer via ``int()``.
            - *buffer* (int): integer number of grid points to keep on each side
              of *crop_extent* (when possible) -- NOT a coordinate distance (cf.
              ``interp_unstructured``'s ``buffer_length``, which is in meters).
              Truncated to an integer via ``int()``.  Counted in *coarsened
              output* points: the operations apply in the order crop -> align ->
              buffer -> coarsen, so the native index window is widened by
              ``buffer * coarsen`` and the strided subsample then keeps
              ``buffer`` of those per edge.  ``buffer=2, coarsen=4`` therefore
              adds 2 points to each edge of the result, not 8.  Fortran
              ``apply_align_buffer_coarsen`` does the same arithmetic.
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

        crop_extent = _resolve_crop_extent(crop_extent,
                                           {'filter_region': filter_region})

        # buffer and coarsen are integer grid-point counts (they feed the index
        # arithmetic below); truncate any float via int() as documented, rather
        # than failing later with an opaque "slice indices must be integers".
        buffer = int(buffer)
        coarsen = int(coarsen)

        if self.unstructured:
            raise NotImplementedError("*** Cannot currently crop unstructured topo")

        if crop_extent is None:
            # only want to coarsen, so this is entire region:
            #crop_extent = [self.x[0],self.x[-1],self.y[0],self.y[-1]]
            crop_extent = self.extent

        xlower,xupper,ylower,yupper = crop_extent

        dx,dy = self.delta
        dx_new = dx*coarsen
        dy_new = dy*coarsen

        # Find crop+coarsen+align index window (shared with the topo_type=4
        # read path so ASCII and NetCDF reads of the same data match exactly).
        idx = _crop_indices(self.x, self.y, crop_extent, coarsen, buffer, align)
        if idx is None:
            # Warned rather than printed so a caller can catch, filter or
            # escalate it; the Fortran reader treats the same condition as
            # fatal (topo_module.f90: "does not overlap topo file", stop 1),
            # so a run that ignores this here will fail there.
            warnings.warn(
                f"crop_extent {list(crop_extent)} does not overlap this "
                f"topography (extent {list(self.extent)}); no crop applied.")
            return None
        ilower, iupper, jlower, jupper = idx

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



def fetch_remote_topo(name_or_url, crop_extent=None, coarsen=1, buffer=0,
                      align=None, nc_params={}, verbose=False):
    r"""Resolve a remote (or local) netCDF DEM into a `Topography`.

    This is the modern one-call "remote DEM -> Topography" path.  It resolves a
    nickname or URL and reads it through the `topo_type=4` reader
    (`Topography.read`, backed by `netcdf_utils.TopoInspector`), so it inherits
    that path's unit handling (a recognized non-meter unit such as `km` is
    converted on read with a warning; a file with no `units` attribute needs
    `assume_units` via `nc_params`), datum handling, fill->NaN conversion, CF
    coordinate/variable detection, and lazy hyperslab windowing.

    :Input:

     - *name_or_url* (str) - a key into `topotools.remote_topo_urls`, or a URL
       (OPeNDAP/THREDDS `dodsC` URLs are read by xarray's netCDF4 backend), or a
       path to a local netCDF file.
     - *crop_extent* ([x1, x2, y1, y2] or None) - requested crop in domain
       coordinates; `None` reads the whole file.  Only the requested hyperslab
       is read from a remote file.
     - *coarsen* (int) - factor to coarsen by (1 = no coarsening).
     - *buffer* (int) - when possible, keep at least this many points outside
       `crop_extent` on each side.
     - *align* ((xalign, yalign) or None) - desired alignment when coarsening;
       see `Topography.crop`.
     - *nc_params* (dict) - options forwarded to the `topo_type=4` reader, e.g.
       `z_var` (elevation variable name) or `assume_units` (unit to assume when
       the file has no `units` attribute).  See `Topography.read`.
     - *verbose* (bool) - if True, print the resolved source.

    :Output:

     - a `topotools.Topography` object.

    Remote-read failures propagate as `OSError`/`RuntimeError` so callers (and
    tests marked `@pytest.mark.remote`) can skip when a server is unavailable.

    Sample usage:

        from clawpack.geoclaw import topotools
        topo = topotools.fetch_remote_topo('etopo22_30sec',
                                           crop_extent=[-126, -122, 46, 49],
                                           coarsen=2, buffer=1, verbose=True)
        topo.write('etopo_sample.tt3', topo_type=3)
    """

    # Resolve a nickname; otherwise treat as a URL or local path.
    if name_or_url in remote_topo_urls:
        url = remote_topo_urls[name_or_url]
    else:
        url = name_or_url

    if verbose:
        print("Will read netCDF data from \n    %s" % url)

    # Set the preprocessing attributes *before* reading: Topography.__init__
    # reads immediately when constructed with a path, which would default these
    # away, so construct empty and read explicitly.
    topo = Topography()
    topo.crop_extent = crop_extent
    topo.coarsen = coarsen
    topo.buffer = buffer
    topo.align = align

    try:
        topo.read(path=url, topo_type=4, nc_params=nc_params)
    except (OSError, RuntimeError):
        # Remote/OPeNDAP servers are flaky; let callers/tests decide to skip.
        raise
    except Exception as e:
        raise RuntimeError(
            "Failed to read remote topo from %s: %s" % (url, e)) from e

    return topo


def read_netcdf(path, zvar=None, extent='all', coarsen=1, return_topo=True,
                return_xarray=False, buffer=0, align=None, verbose=False):

    r"""Deprecated: read a netCDF DEM into a Topography and/or xarray.Dataset.

    .. deprecated::
        Use :func:`fetch_remote_topo` (or ``Topography.read(topo_type=4)``)
        instead.  This is now a thin wrapper over :func:`fetch_remote_topo`; the
        standalone ``netCDF4``-based reader it used to contain has been removed
        in favor of the modern ``topo_type=4`` read path (unit checking, datum,
        fill->NaN, CF coordinate/variable detection, lazy hyperslab windowing).

    The legacy signature is preserved:

     - *path* (str) - nickname (key of ``remote_topo_urls``), URL, or local file.
     - *zvar* (str) - elevation variable name; mapped to ``nc_params['z_var']``.
     - *extent* - ``[x1,x2,y1,y2]`` requested crop, or ``'all'`` for whole file.
     - *coarsen* (int) - coarsening factor (1 = none).
     - *return_topo* (bool) - if True, include a ``Topography`` in the result.
     - *return_xarray* (bool) - if True, include an ``xarray.Dataset``.
     - *buffer* (int) - points to keep outside the crop on each side.
     - *align* (tuple) - alignment when coarsening; see ``Topography.crop``.

    :Output:
     - a ``Topography``, an ``xarray.Dataset``, or a ``(topo, ds)`` tuple,
       depending on ``return_topo`` / ``return_xarray`` (unchanged contract).
    """

    warnings.warn(
        "topotools.read_netcdf is deprecated; use "
        "topotools.fetch_remote_topo (or Topography.read(topo_type=4)) instead.",
        DeprecationWarning, stacklevel=2)

    assert (type(coarsen) is int) and (coarsen >= 1), \
        '*** coarsen must be a positive integer'

    # Map the legacy arguments onto the modern helper.
    crop_extent = None if (isinstance(extent, str) and extent == 'all') \
        else extent
    nc_params = {}
    if zvar is not None:
        nc_params['z_var'] = zvar

    topo = fetch_remote_topo(path, crop_extent=crop_extent, coarsen=coarsen,
                             buffer=buffer, align=align, nc_params=nc_params,
                             verbose=verbose)

    output = None
    if return_topo:
        output = topo

    if return_xarray:
        import xarray
        # Rebuild an xarray.Dataset from the resulting Topography so the legacy
        # return contract is unchanged.  Z has shape (len(y), len(x)).
        xarray_ds = xarray.Dataset({'z': (('lat', 'lon'), topo.Z)},
                                   coords={'lon': topo.x, 'lat': topo.y})
        if output is None:
            output = xarray_ds
        else:
            output = (topo, xarray_ds)

    return output
