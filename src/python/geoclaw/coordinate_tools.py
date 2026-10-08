#!/usr/bin/env python
# encoding: utf-8
r"""
Format- and product-neutral coordinate helpers for GeoClaw input.

This module holds the longitude/geometry primitives that are shared by *all*
input products (topography, moving topography / dtopo, and gridded met forcing)
and *all* file formats (ASCII types 1/2/3 and NetCDF type 4).  Keeping them here
-- rather than in ``netcdf_utils`` -- gives every path a single implementation
to import, so the antimeridian-wrap logic is never copied per format.

Contents
--------
is_geographic_lon
    Decide whether an x axis is a geographic longitude (so the +/-360 wrap is
    meaningful) or a projected / Cartesian axis (which never wraps).
_compute_lon_entries
    Pure geometry: given a file's longitude span and a requested domain span,
    return the ``(file_crop_min, file_crop_max, lon_offset)`` tuples needed to
    cover the domain -- one tuple normally, two when the request straddles the
    file's longitude cut (dateline).
"""

from __future__ import annotations

from typing import Optional

import numpy


# Length unit strings (lower-cased) on an x axis that mark it as a projected /
# rectilinear (non-geographic) grid, for which the 0-360 longitude wrap is
# meaningless and must be skipped.  Detection defaults to geographic; only
# positive evidence (these units, a projection standard_name, or values
# outside the plausible degree range) disables the wrap.
_PROJECTED_LENGTH_UNITS: frozenset[str] = frozenset({
    'm', 'meter', 'meters', 'metre', 'metres',
    'km', 'kilometer', 'kilometers', 'kilometre', 'kilometres',
})
_PROJECTED_STANDARD_NAMES: frozenset[str] = frozenset({
    'projection_x_coordinate', 'projection_y_coordinate',
})

# Positive evidence (lower-cased) that an x axis IS geographic longitude: CF
# degree units or a longitude standard_name.  Used to distinguish a file that is
# *definitely* geographic from one that merely lacks metadata (the ambiguous
# default) -- the coordinate-system gate only hard-errors on positive evidence.
_GEOGRAPHIC_LON_UNITS: frozenset[str] = frozenset({
    'degrees_east', 'degree_east', 'degrees_e', 'degree_e', 'degreese',
    'degrees', 'degree', 'degrees_north', 'degree_north',  # generic degree axes
})
_GEOGRAPHIC_STANDARD_NAMES: frozenset[str] = frozenset({
    'longitude', 'grid_longitude',
})


def classify_lon_axis(
    units: Optional[str] = None,
    std_name: Optional[str] = None,
    lon_min: Optional[float] = None,
    lon_max: Optional[float] = None,
) -> str:
    """Classify an x axis as ``'projected'``, ``'geographic'``, or ``'unknown'``.

    Only *positive* evidence promotes an axis out of ``'unknown'``:
      - ``'projected'`` -- a length unit (``m``/``km``), a projection
        ``standard_name``, or coordinate values outside the plausible degree
        range ``[-360, 360]``.
      - ``'geographic'`` -- a CF degree unit or a longitude ``standard_name``.

    An axis with no such metadata (e.g. an ASCII file, or a small unit-less
    Cartesian grid) is ``'unknown'``.  This three-way result lets the
    coordinate-system gate (:func:`resolve_wrap`) hard-error only on a genuine,
    positively-evidenced conflict and never on an ambiguous default.
    """
    u = str(units or '').strip().lower()
    s = str(std_name or '').strip().lower()
    if (u in _PROJECTED_LENGTH_UNITS
            or s in _PROJECTED_STANDARD_NAMES
            or (lon_min is not None and lon_min < -360.0 - 1e-6)
            or (lon_max is not None and lon_max > 360.0 + 1e-6)):
        return 'projected'
    if u in _GEOGRAPHIC_LON_UNITS or s in _GEOGRAPHIC_STANDARD_NAMES:
        return 'geographic'
    return 'unknown'


def is_geographic_lon(
    lon_min: float,
    lon_max: float,
    units: Optional[str] = None,
    std_name: Optional[str] = None,
) -> bool:
    """Return True if an x axis should be treated as geographic longitude.

    The axis is assumed geographic unless there is positive evidence it is a
    projected / Cartesian (e.g. beta-plane) grid: a length unit (``m``/``km``),
    a projection ``standard_name``, or coordinate values outside the plausible
    degree range ``[-360, 360]``.  Antimeridian wrapping (+/-360) is only ever
    applied when this returns True; a Cartesian grid is left untouched.

    ``units``/``std_name`` are optional (ASCII files carry no such metadata);
    when omitted, only the coordinate-magnitude test applies.
    """
    return classify_lon_axis(units, std_name, lon_min, lon_max) != 'projected'


def resolve_wrap(coordinate_system: Optional[int], nature: str) -> bool:
    """Decide whether antimeridian wrapping is allowed, gated on the run's
    ``coordinate_system`` -- the single authority (GeoClaw ``geodata``:
    ``1`` = Cartesian/meters, ``2`` = lon-lat sphere).

    *nature* is the file x-axis classification from :func:`classify_lon_axis`.

    Returns ``allow_wrap``.  Raises ``ValueError`` on a genuine mismatch (a
    positively geographic file under a Cartesian run, or a positively projected
    file under a geographic run) -- symmetric, per the design decision that
    mixing coordinate systems is always a setup error.

    ``coordinate_system is None`` means the run's system is unknown at this call
    site (e.g. the Python early-check was not wired a value); wrapping then
    falls back to the per-file heuristic (geographic unless positively
    projected), preserving legacy behavior with no cross-check.
    """
    if coordinate_system == 1:      # Cartesian / UTM: never wrap.
        if nature == 'geographic':
            raise ValueError(
                "Coordinate-system mismatch: a geographic (lon/lat) input file "
                "was supplied for a Cartesian run (geodata.coordinate_system == "
                "1). Longitude wrapping is disabled for Cartesian coordinates; "
                "use a projected input file or set coordinate_system = 2.")
        return False
    if coordinate_system == 2:      # lon-lat sphere: wrapping is meaningful.
        if nature == 'projected':
            raise ValueError(
                "Coordinate-system mismatch: a projected/Cartesian input file "
                "was supplied for a geographic run (geodata.coordinate_system == "
                "2). Use a lon/lat input file or set coordinate_system = 1.")
        return True
    # Unknown coordinate_system: legacy per-file behavior, no cross-check.
    return nature != 'projected'


def _compute_lon_entries(
    file_lon_min: float,
    file_lon_max: float,
    domain_lon_min: float,
    domain_lon_max: float,
    max_gap: float = 1e-10,
    allow_wrap: bool = True,
) -> list[tuple[float, float, float]]:
    """
    Compute (file_crop_min, file_crop_max, lon_offset) tuples needed to cover
    [domain_lon_min, domain_lon_max] from a file with lons in
    [file_lon_min, file_lon_max].

    lon_offset is the scalar added to file coordinates to produce domain
    coordinates: x_domain = x_file + lon_offset.  For NetCDF this is the
    descriptor ``lon_wrap_offset``; for ASCII it is folded into ``x_shift``.

    Returns 1 tuple if a single offset suffices, 2 tuples if the domain
    straddles the file's cut point.

    max_gap controls how much under-coverage is tolerated.  The default
    (1e-10) is tight enough to catch genuine gaps.  Pass the file's grid
    spacing to allow for the half-cell gap at the dateline that near-global
    files (e.g. GEBCO) have between their last and first longitude columns.

    allow_wrap enables the +/-360 wrap candidate offsets (the default, for a
    geographic longitude axis).  Pass False for a non-geographic x axis
    (projected meters, etc.), which never wraps: only the identity offset 0
    is considered.

    Raises ValueError if the file cannot cover the requested domain even
    with all candidate offsets, or if the remaining uncovered gap exceeds
    max_gap.
    """
    candidate_offsets = [0.0, 360.0, -360.0] if allow_wrap else [0.0]
    entries: list[tuple[float, float, float]] = []
    total_coverage = 0.0

    for offset in candidate_offsets:
        shifted_min = file_lon_min + offset
        shifted_max = file_lon_max + offset
        intersect_min = max(domain_lon_min, shifted_min)
        intersect_max = min(domain_lon_max, shifted_max)
        width = intersect_max - intersect_min
        if width > 1e-10:
            file_crop_min = intersect_min - offset
            file_crop_max = intersect_max - offset
            # Clamp to file extent (guards against floating-point overshoot)
            file_crop_min = max(file_crop_min, file_lon_min)
            file_crop_max = min(file_crop_max, file_lon_max)
            entries.append((file_crop_min, file_crop_max, offset))
            total_coverage += width

    if not entries:
        raise ValueError(
            f"File longitude range [{file_lon_min}, {file_lon_max}] cannot cover "
            f"domain [{domain_lon_min}, {domain_lon_max}] with candidate offsets "
            f"{candidate_offsets}."
        )

    # One-sided check: overcoverage is harmless; under-coverage beyond max_gap
    # indicates that the file genuinely cannot cover the requested domain.
    gap = (domain_lon_max - domain_lon_min) - total_coverage
    if gap > max_gap:
        raise ValueError(
            f"File longitudes [{file_lon_min}, {file_lon_max}] cannot "
            f"cover requested domain [{domain_lon_min}, {domain_lon_max}]. "
            f"Gap of {gap:.6f} degrees exceeds tolerance {max_gap:.6f}."
        )

    return entries


class EmptyCropWindow(ValueError):
    r"""A crop request overlaps the coordinate range but contains no grid point.

    Distinct from a genuine non-overlap (which :func:`crop_indices` reports by
    returning ``None``, and whose callers have a documented full-file
    fall-back): a window narrower than one cell, falling strictly between two
    grid points, is never what the user meant, so it is raised rather than
    quietly widened.  Subclasses ``ValueError`` so existing ``except
    ValueError`` handling is unaffected; callers that can add axis names or
    grid spacings to the message catch it and re-raise.
    """


def crop_indices(coords, lo, hi, delta, coarsen=1, buffer=0, align=None):
    r"""Half-open slice bounds ``(lower, upper)`` for a 1-D crop + coarsen.

    ``coords[lower:upper:coarsen]`` is the cropped, coarsened, buffered,
    align-adjusted lattice covering the inclusive request ``[lo, hi]``.  This is
    the single index-window implementation shared by the topography and moving-
    topography (dtopo) crop paths, for both ASCII and NetCDF products.

    :Input:
     - *coords* (ndarray) - 1-D ascending coordinate array with uniform spacing.
     - *lo*, *hi* (float) - inclusive requested crop bounds (domain coords).
     - *delta* (float) - grid spacing of *coords* (``coords[1] - coords[0]``).
     - *coarsen* (int) - coarsening (subsampling) factor; 1 = none.
     - *buffer* (int) - integer number of grid points to keep on each side (NOT
       a coordinate distance).  Scaled by *coarsen* so it counts points on the
       original, pre-coarsen lattice.
     - *align* (float or None) - alignment target for coarsening; when given and
       ``coarsen > 1`` the low index is shifted so ``(coords[lower] - align) /
       (delta*coarsen)`` is as close to an integer as the lattice allows.

    Returns ``(lower, upper)`` (Python ints, half-open) or ``None`` when
    ``[lo, hi]`` does not overlap *coords* -- callers preserve their own
    "does not overlap" handling.

    :Raises:
     - :class:`EmptyCropWindow` when ``[lo, hi]`` lies inside the coordinate
       range but between two grid points, so the window is empty.
    """
    coords = numpy.asarray(coords)
    n = len(coords)
    try:
        lower = (coords >= lo).nonzero()[0][0]
        upper = (coords <= hi).nonzero()[0][-1]
    except IndexError:
        return None
    if upper < lower:
        raise EmptyCropWindow(
            f"crop bounds [{lo}, {hi}] lie between grid points and contain no "
            f"data: the grid spacing is {delta}. Widen the crop to at least "
            f"one cell, or use buffer= to include the surrounding points.")

    delta_new = delta * coarsen

    # Shift the low index for alignment when coarsening (mirrors the original
    # Topography.crop math exactly: search the first `coarsen` points for the
    # best-aligned start, then trim the high index to a whole coarsen stride).
    if (coarsen > 1) and (align is not None):
        vs = numpy.array([coords[lower + i] for i in range(coarsen)])
        offsets = (vs - align) / delta_new
        offsets_frac = offsets - numpy.round(offsets)
        ioffset = numpy.argmin(abs(offsets_frac))
        lower = lower + ioffset
        upper = upper - numpy.remainder(upper - lower, coarsen)

    # Buffer (integer grid-point count), clamped to the array limits.
    lower = numpy.maximum(0, lower - buffer * coarsen)
    upper = numpy.minimum(n - 1, upper + buffer * coarsen) + 1

    return int(lower), int(upper)


def axis_file_slice(coord_full, descending, lo, hi, step, n):
    r"""Map an ascending window ``[lo:hi:step]`` to a positive-stride slice.

    Used by every NetCDF input reader: :func:`crop_indices` returns bounds into
    an *ascending* view of a file axis, but NetCDF/xarray lazy indexing requires
    a **positive** step into the *file-order* axis.  For an axis stored
    descending (e.g. latitude N->S), the ascending sample indices
    ``lo, lo+step, ..., lo+(m-1)*step`` map to file indices ``n-1-(that)``, whose
    minimum is ``f0``; reading ``slice(f0, f0+m*step, step)`` returns those same
    ``m`` samples in file (descending) order, to be flipped to ascending in
    memory afterward.

    Returns ``(file_slice, coord_subset, flip)`` where
    ``coord_subset == coord_full[file_slice]`` and ``flip`` is True when the
    subset (and the corresponding data axis) must be reversed to be ascending.
    """
    if not descending:
        sl = slice(lo, hi, step)
        return sl, coord_full[sl], False
    m = len(range(lo, hi, step))
    f0 = n - 1 - (lo + (m - 1) * step)
    sl = slice(f0, f0 + m * step, step)
    return sl, coord_full[sl], True
