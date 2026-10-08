#!/usr/bin/env python
# encoding: utf-8
r"""
Extensible format-reader seam for GeoClaw gridded input.

Every gridded input product (topography and, in later work, moving topography /
dtopo and gridded met forcing) needs to ingest a lat/lon-or-Cartesian grid from
one of several on-disk formats.  Rather than let each product hard-code a
``topo_type``/format ``if/elif`` ladder, the format-specific *reading* lives in
a small family of :class:`GriddedReader` adapters registered here.  A product
asks :func:`get_reader` for the reader that handles a given ``topo_type`` and
calls :meth:`GriddedReader.read_window`; all shared post-processing (unit /
fill / preprocessing handling) stays in the product.

This is the single place a new format (GeoTIFF today, zarr later) is added: write
an adapter and register it.  Formats the Fortran core cannot read natively
(GeoTIFF, zarr) are *ingest-only* -- read into memory here, then written back out
as NetCDF (type 4) or ASCII (type 3) for the simulation.

Adapters are intentionally thin wrappers over the existing, tested read paths;
this module concentrates the *dispatch*, not new I/O logic.
"""

from __future__ import annotations

import dataclasses
from typing import Optional

import numpy

from clawpack.geoclaw import coordinate_tools


# Sentinel distinguishing "reader did not produce this field" from a produced
# ``None`` (e.g. the NetCDF reader always resolves a datum, which may be None).
_UNSET = object()


@dataclasses.dataclass
class GridData:
    r"""Raw grid returned by a :class:`GriddedReader`.

    ``x``/``y`` are 1-D coordinate arrays (ascending, S->N for ``y``); ``Z`` is
    the 2-D field in ``(y, x)`` order.  ``delta`` and ``datum`` are optional
    extras a reader may resolve; left unset otherwise so the caller keeps its
    existing value.
    """
    x: numpy.ndarray
    y: numpy.ndarray
    Z: numpy.ndarray
    delta: Optional[tuple] = None
    datum: object = _UNSET


class GriddedReader:
    r"""Base class for a format adapter.

    Subclasses declare the ``topo_types`` they handle and implement
    :meth:`read_window`, which returns a :class:`GridData` for the *topo*
    context (its ``path``, ``topo_type`` and preprocessing attributes such as
    ``crop_extent``/``coarsen``/``buffer``/``align``).
    """

    #: ``abs(topo_type)`` values this reader handles.
    topo_types: frozenset = frozenset()

    #: True when :meth:`read_window` already applies the caller's
    #: crop/coarsen/buffer/align while reading, so the caller must NOT run its
    #: own crop afterward (doing so would coarsen twice).
    applies_preprocessing: bool = False

    def read_window(self, topo, *, nc_params=None) -> GridData:
        raise NotImplementedError


class ASCIIReader(GriddedReader):
    r"""ASCII DEM formats: type 1 (x y z columns) and types 2/3 (headered grid).

    Reads the whole file; the caller applies crop/coarsen/buffer/align to the
    in-memory arrays afterward (``applies_preprocessing`` is False).
    """

    topo_types = frozenset({1, 2, 3})

    def read_window(self, topo, *, nc_params=None) -> GridData:
        if abs(topo.topo_type) == 1:
            return self._read_type1(topo)
        return self._read_type23(topo)

    def _read_type1(self, topo) -> GridData:
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
            topo.crop_extent is not None
            or topo.coarsen != 1
            or topo.buffer != 0
            or topo.align is not None
            or topo.x_shift != 0.0
            or topo.y_shift != 0.0
            or topo.z_shift != 0.0
            or topo.negate_z
        )
        if _preprocessing_requested:
            raise NotImplementedError(
                "Preprocessing attributes (crop_extent, coarsen, buffer, align, "
                "shift) are not supported for topo_type=1. Convert to type 2/3/4 first:\n"
                "  topo_raw = Topography(path=self.path, topo_type=1)\n"
                "  topo_raw.read()\n"
                "  topo_raw.write('converted.tt2', topo_type=2)"
            )
        data = numpy.loadtxt(topo.path)
        N = [0, 0]
        y0 = data[0, 1]
        for (n, y) in enumerate(data[1:, 1]):
            if y != y0:
                N[1] = n + 1
                break
        N[0] = data.shape[0] // N[1]

        x = data[:N[1], 0]
        y = data[::N[1], 1]
        Z = numpy.flipud(data[:, 2].reshape(N))
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        return GridData(x=x, y=y, Z=Z, delta=(dx, dy))

    def _read_type23(self, topo) -> GridData:
        # read_header sets topo._extent, _x, _y, _delta, grid_registration.
        N = topo.read_header()
        if abs(topo.topo_type) == 2:
            # Single column, reshape it.
            Z = numpy.loadtxt(topo.path, skiprows=6).reshape(N[1], N[0])
            Z = numpy.flipud(Z)
        else:  # abs(topo.topo_type) == 3
            # Read starting at the top-left corner.
            Z = numpy.flipud(numpy.loadtxt(topo.path, skiprows=6))
        return GridData(x=topo._x, y=topo._y, Z=Z, delta=topo._delta)


class NetCDFReader(GriddedReader):
    r"""CF-compliant NetCDF DEM (type 4), read via netcdf_utils.TopoInspector.

    Unlike :class:`ASCIIReader` this applies the caller's crop/coarsen/buffer/
    align *while reading*, pushing the window down to the NetCDF hyperslab so a
    global DEM is never materialized in full (``applies_preprocessing`` is
    True).  The window comes from the same ``_crop_indices`` routine
    ``Topography.crop`` uses, so an ASCII read + crop and a NetCDF read of the
    same data produce identical grids.
    """

    topo_types = frozenset({4})
    applies_preprocessing = True

    def read_window(self, topo, *, nc_params=None) -> GridData:
        if nc_params is None:
            nc_params = {}
        import warnings

        from clawpack.geoclaw import netcdf_utils as _ncutils
        # Imported as a module (not `from ... import _crop_indices`) so the
        # lookup happens at call time: the crop-pushdown tripwire test
        # monkeypatches topotools._crop_indices and must see the substitution.
        from clawpack.geoclaw import topotools as _topotools
        from clawpack.geoclaw.topotools import extract_datum
        from clawpack.geoclaw.units import (
            convert as _units_convert,
            GEOCLAW_NETCDF_UNITS as _NC_UNITS,
        )

        # Explicit variable name override via nc_params (back-compat 'z_var').
        _var_hint = nc_params.get('z_var', None)
        # Opt-in escape hatch for a file with no 'units'; units are otherwise
        # required, never silently assumed.
        _assume_units = nc_params.get('assume_units', None)
        # Opt-out of the post-conversion magnitude sanity check.
        _skip_sanity = nc_params.get('skip_sanity_check', False)

        datum = _UNSET
        with _ncutils.TopoInspector(
            topo.path, var_name=_var_hint, assume_units=_assume_units,
            skip_sanity_check=_skip_sanity,
        ) as inspector:
            if inspector.var_name is None:
                inspector.var_name = inspector._find_topo_var_name()

            # CF coordinate metadata (no fill-value data scan here; fill values
            # are handled by the caller via NaN replacement).
            _meta = inspector.inspect(inspector.var_name)
            # Check/record units; warns if conversion is needed.
            _source_units = inspector._check_topo_units()

            ds = inspector.ds
            _x_name = _meta.x_name
            _y_name = _meta.y_name
            _var_name = inspector.var_name

            # Optional vertical datum (informational only; no transformation).
            datum = extract_datum(ds[_var_name].attrs, ds.attrs)

            # Load 1-D coordinate arrays in full (O(nx)+O(ny), cheap); any N->S
            # flip is applied below after windowing.
            _lon_full = numpy.asarray(ds[_x_name].values, dtype=float)
            _lat_full = numpy.asarray(ds[_y_name].values, dtype=float)

            # Load variable, squeezing singleton non-spatial dims.
            _da = ds[_var_name]
            for _dim in list(_da.dims):
                if _dim not in (_x_name, _y_name):
                    if _da.sizes[_dim] == 1:
                        _da = _da.isel({_dim: 0})
                    else:
                        raise ValueError(
                            f"NetCDF variable '{_var_name}' has non-singleton "
                            f"dimension '{_dim}' (size {_da.sizes[_dim]}).  "
                            f"Cannot load as static topography."
                        )

            # Transpose to (lat, lon) = (y, x) order expected by Topography.
            _da = _da.transpose(_y_name, _x_name)

            # Push the crop+coarsen+align window down to xarray's lazy indexing
            # so the NetCDF backend reads ONLY the final hyperslab.  Otherwise
            # `_da.values` materializes the whole variable (a global DEM is many
            # GB, and CF fill decoding promotes it to float), which is
            # prohibitively slow and can exhaust memory even when the caller
            # asked for a small, coarsened subset.  The window is computed on
            # the cheap 1-D coordinate arrays via _crop_indices -- the SAME
            # routine crop() uses -- so this read matches an ASCII read + crop()
            # of the same data exactly; the caller then skips its post-read
            # crop() (applies_preprocessing is True).
            _nx = _lon_full.size
            _ny = _lat_full.size
            _c = max(int(topo.coarsen), 1)

            # _crop_indices requires ascending coords; build ascending views of
            # the file axes (lon is usually ascending, lat is often N->S).
            _lon_desc = _lon_full[0] > _lon_full[-1]
            _lat_desc = not _meta.y_increasing
            _lon_asc = _lon_full[::-1] if _lon_desc else _lon_full
            _lat_asc = _lat_full[::-1] if _lat_desc else _lat_full

            if topo.crop_extent is not None:
                _ce = list(topo.crop_extent)
            else:
                # whole file (mirrors crop()'s crop_extent=self.extent)
                _ce = [_lon_asc[0], _lon_asc[-1], _lat_asc[0], _lat_asc[-1]]

            _idx = _topotools._crop_indices(_lon_asc, _lat_asc, _ce, _c,
                                            int(topo.buffer), topo.align)
            if _idx is None:
                # crop_extent misses the file: fall back to the full grid at
                # native resolution (mirrors crop() returning None -> no-op),
                # rather than coarsening the whole file.
                warnings.warn(
                    f"crop_extent {list(_ce)} did not overlap {topo.path}; the "
                    f"full file was read uncropped at native resolution. The "
                    f"Fortran reader treats a non-overlapping crop as fatal.")
                _il, _iu, _jl, _ju = 0, _nx, 0, _ny
                _step = 1
            else:
                _il, _iu, _jl, _ju = _idx
                _step = _c

            # Map each ascending [lo:hi:step] window to a positive-stride slice
            # into the FILE-order axis; descending axes are flipped in memory
            # afterward.
            _x_slice, _lon_vals, _flip_x = coordinate_tools.axis_file_slice(
                _lon_full, _lon_desc, _il, _iu, _step, _nx)
            _y_slice, _lat_vals, _flip_y = coordinate_tools.axis_file_slice(
                _lat_full, _lat_desc, _jl, _ju, _step, _ny)

            _da = _da.isel({_y_name: _y_slice, _x_name: _x_slice})
            _z_vals = numpy.asarray(_da.values, dtype=float)

            # Flip descending axes to ascending (S->N, W->E) in memory.
            if _flip_y:
                _lat_vals = _lat_vals[::-1]
                _z_vals = _z_vals[::-1, :]
            if _flip_x:
                _lon_vals = _lon_vals[::-1]
                _z_vals = _z_vals[:, ::-1]

            # Unit conversion if source is not already meters.
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
                    var_name=_var_name, path=str(topo.path),
                )

            # Decoded fill values are already NaN (xarray mask_and_scale); NaN
            # is the in-memory missing-data representation, so they pass through.

        return GridData(x=_lon_vals, y=_lat_vals, Z=_z_vals, datum=datum)


class GeoTIFFReader(GriddedReader):
    r"""GeoTIFF DEM (type 5) via GDAL.  Ingest-only: transcode to type 3/4 for Fortran."""

    topo_types = frozenset({5})

    def read_window(self, topo, *, nc_params=None) -> GridData:
        try:
            import gdal
        except ImportError as e:
            print("Reading GeoTIFF files requires GDAL.")
            raise e

        data = gdal.Open(topo.path)
        z = data.GetRasterBand(1).ReadAsArray()
        transform = data.GetGeoTransform()
        x_origin = transform[0]
        y_origin = transform[3]
        dx = transform[1]
        dy = -transform[5]

        Z = numpy.flipud(z)
        x = numpy.linspace(x_origin,
                           x_origin + (z.shape[0] - 1) * dx, z.shape[0])
        y = numpy.linspace(y_origin - (z.shape[1] - 1) * dy,
                           y_origin, z.shape[1])
        return GridData(x=x, y=y, Z=Z)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_REGISTRY: dict[int, GriddedReader] = {}


def register_reader(reader: GriddedReader) -> None:
    r"""Register *reader* for each ``topo_type`` in ``reader.topo_types``."""
    for topo_type in reader.topo_types:
        _REGISTRY[int(topo_type)] = reader


def get_reader(topo_type) -> GriddedReader:
    r"""Return the registered :class:`GriddedReader` for *topo_type*.

    Dispatch is by ``abs(topo_type)`` (the sign only flips the elevation
    convention, handled by the caller).  Raises ``IOError`` for an unknown type.
    """
    try:
        return _REGISTRY[abs(int(topo_type))]
    except (KeyError, TypeError, ValueError):
        raise IOError("Unrecognized topo_type: %s" % topo_type)


register_reader(ASCIIReader())
register_reader(NetCDFReader())
register_reader(GeoTIFFReader())
