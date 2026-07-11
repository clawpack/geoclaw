#!/usr/bin/env python
# encoding: utf-8
r"""
NetCDF utilities for GeoClaw NetCDF input.

Provides inspector classes that inspect NetCDF files for coordinate
metadata, unit consistency, and fill values without loading data arrays.
The inspectors produce metadata dataclasses consumed by DescriptorWriter
(not yet implemented — pending confirmation of Fortran topo.data parser).

Classes
-------
NetCDFInspector
    Base class: opens a file, discovers coordinate variables, detects lon/lat
    conventions, resolves fill value, validates crop bounds.

TopoInspector(NetCDFInspector)
    Adds bathymetry-specific checks: fill values within crop region (fatal),
    unit verification against GEOCLAW_NETCDF_UNITS['topo'].

MetInspector(NetCDFInspector)
    Adds met-forcing checks: multi-variable grid/time consistency, unit
    verification and conversion to contract units, CF time -> seconds offset,
    ensemble dimension detection.

CFNormalizer
    Standalone utility: renames coordinates to CF standard names, adds/fixes
    standard_name / axis / units / _FillValue attributes, resolves
    _FillValue vs missing_value conflicts. Does not resample or reproject.

Metadata dataclasses
--------------------
FileMetadata, TopoMetadata, MetVariableInfo, MetMetadata
"""

from __future__ import annotations

import contextlib
import dataclasses
import importlib.util
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import xarray as xr

from clawpack.geoclaw.units import GEOCLAW_NETCDF_UNITS, convert as units_convert


# ---------------------------------------------------------------------------
# CF unit alias table
# Maps CF-style unit strings to the abbreviations understood by units.convert()
# ---------------------------------------------------------------------------

# For each contract role, the set of CF unit strings that mean "already correct"
_CONTRACT_UNIT_CF_ALIASES: dict[str, frozenset[str]] = {
    'm':   frozenset({'m', 'meter', 'meters', 'metre', 'metres'}),
    'm/s': frozenset({'m/s', 'm s-1', 'meters per second', 'metres per second',
                      'meter per second', 'metre per second', 'm s**-1'}),
    'Pa':  frozenset({'Pa', 'pa', 'pascal', 'pascals', 'Pascal', 'Pascals',
                      'N/m2', 'N m-2', 'kg m-1 s-2'}),
    's':   frozenset({'s', 'sec', 'second', 'seconds'}),
}

# Maps CF unit strings to units.py abbreviations for conversion
_CF_TO_UNITS_PY: dict[str, str] = {
    # length
    'm': 'm', 'meter': 'm', 'meters': 'm', 'metre': 'm', 'metres': 'm',
    'cm': 'cm', 'centimeter': 'cm', 'centimeters': 'cm',
    'km': 'km', 'kilometer': 'km', 'kilometers': 'km',
    # speed
    'm/s': 'm/s', 'm s-1': 'm/s', 'm s**-1': 'm/s',
    'meters per second': 'm/s', 'metres per second': 'm/s',
    'knots': 'knots', 'kt': 'knots', 'kts': 'knots',
    'km/h': 'km/h', 'km h-1': 'km/h',
    'mph': 'mph',
    # pressure
    'Pa': 'Pa', 'pa': 'Pa', 'pascal': 'Pa', 'pascals': 'Pa',
    'hPa': 'hPa', 'hectopascal': 'hPa', 'hectopascals': 'hPa',
    'mbar': 'mbar', 'mb': 'mbar', 'millibar': 'mbar', 'millibars': 'mbar',
    # time
    's': 's', 'sec': 's', 'second': 's', 'seconds': 's',
    'min': 'min', 'minute': 'min', 'minutes': 'min',
    'h': 'h', 'hr': 'h', 'hour': 'h', 'hours': 'h',
    'days': 'days', 'day': 'days',
}


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


def _normalize_cf_unit(cf_unit: str) -> Optional[str]:
    """Return the units.py abbreviation for *cf_unit*, or None if unknown."""
    return _CF_TO_UNITS_PY.get(cf_unit.strip())


def _unit_matches_contract(cf_unit: str, contract_unit: str) -> bool:
    """Return True if *cf_unit* is equivalent to *contract_unit*."""
    aliases = _CONTRACT_UNIT_CF_ALIASES.get(contract_unit, frozenset())
    return cf_unit.strip() in aliases


# units.py time-unit abbreviations (the canonical forms _CF_TO_UNITS_PY maps
# time strings onto): used to validate that a bare time-axis units attribute
# really is a duration before scaling it to seconds.
_TIME_UNIT_ABBREVS: frozenset[str] = frozenset({'s', 'min', 'h', 'days'})


def _cf_time_units_to_seconds_factor(cf_unit: str) -> float:
    """Multiplicative factor converting a value in *cf_unit* to seconds.

    *cf_unit* must be a bare CF duration unit (``seconds``, ``minutes``,
    ``hours``, ``days`` and their aliases).  Raises ValueError if it is not a
    recognised time unit -- callers must never silently assume seconds for an
    unrecognised or non-time unit string.
    """
    canonical = _normalize_cf_unit(cf_unit)
    if canonical is None or canonical not in _TIME_UNIT_ABBREVS:
        raise ValueError(
            f"Unrecognised time units {cf_unit!r}; expected a CF duration unit "
            f"such as 'seconds', 'minutes', 'hours', or 'days'."
        )
    return float(units_convert(1.0, canonical, 's'))


@contextlib.contextmanager
def suppress_netcdf4_shape_warning():
    r"""Silence netCDF4-python's spurious NumPy 2.5 shape-assignment warning.

    ``netCDF4.Variable.__setitem__`` reshapes the incoming array in place via
    ``data.shape = ...`` (see ``netCDF4/_netCDF4.pyx``), which NumPy >= 2.5
    now deprecates in favor of ``np.reshape``. This is entirely internal to
    netCDF4-python (as of 1.7.4, the latest release) and not something callers
    can avoid, so wrap ``var[:] = ...`` assignments in this context manager
    rather than letting the warning leak into test output. Safe to remove
    once netCDF4-python ships a fix upstream.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Setting the shape on a NumPy array has been deprecated",
            category=DeprecationWarning,
        )
        yield


# ---------------------------------------------------------------------------
# Metadata dataclasses
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class FileMetadata:
    """Coordinate and convention metadata for a single NetCDF file."""
    source_file: Path
    x_name: str
    y_name: str
    time_name: Optional[str]
    lon_wrap: Optional[int]      # 180 -> [-180,180]; 360 -> [0,360]; None -> non-geographic
    y_increasing: bool           # True if y increases with array index
    dim_order: list[str]         # canonical role names, e.g. ['y', 'x']
    fill_value: Optional[float]  # resolved from _FillValue / missing_value
    crop_bounds: Optional[tuple[float, float, float, float]]  # x0,x1,y0,y1


@dataclasses.dataclass
class TopoMetadata(FileMetadata):
    """FileMetadata plus topo-specific fields."""
    var_name: str
    source_units: str    # units as found in file
    fill_action: str     # 'abort' (only value for topo; fill = fatal)
    lon_wrap_offset: float = 0.0  # scalar Fortran adds to file coords: x_domain = x_file + lon_wrap_offset


@dataclasses.dataclass
class DTopoMetadata(FileMetadata):
    """FileMetadata plus dtopo-specific fields.

    The time axis is collapsed to (t0, dt) in simulation seconds: the
    Fortran dtopo machinery requires uniform time spacing and reconstructs
    times as t0 + k*dt, so Fortran never parses the file's time variable.
    """
    var_name: str        # the dZ (deformation) variable
    t0: float            # first time, simulation seconds
    dt: float            # uniform time step, seconds (0.0 when mt == 1)
    mt: int              # number of times (informational; Fortran uses dims)
    lon_wrap_offset: float = 0.0  # scalar Fortran adds: x_domain = x_file + offset


@dataclasses.dataclass
class MetVariableInfo:
    """Maps one NetCDF variable to its GeoClaw role."""
    var_name: str
    geoclaw_role: str    # 'wind_u', 'wind_v', 'pressure'
    source_units: str    # units as found in file


@dataclasses.dataclass
class MetMetadata(FileMetadata):
    """FileMetadata plus met-forcing-specific fields."""
    variables: list[MetVariableInfo]
    time_offset: float   # seconds; Fortran adds this to times read from file
    fill_action: str     # 'abort' or 'warn'


# ---------------------------------------------------------------------------
# Base inspector
# ---------------------------------------------------------------------------

class NetCDFInspector:
    """
    Open a NetCDF file and inspect its coordinate metadata.

    No data arrays are loaded; only coordinate values and variable attributes
    are accessed.  Dask-lazy chunking is enabled so that any accidental
    downstream `.compute()` calls are bounded.

    Parameters
    ----------
    path : str or Path
        Path to the NetCDF file.
    crop_bounds : tuple of four floats, optional
        (lon_min, lon_max, lat_min, lat_max) in the same convention as the
        file.  When provided, bounds are validated against file extent.
    """

    def __init__(
        self,
        path: str | Path,
        crop_bounds: Optional[tuple[float, float, float, float]] = None,
    ) -> None:
        self.path = Path(path)
        self.crop_bounds = crop_bounds
        # Activate Dask-lazy chunking if dask is available; fall back to
        # netCDF4 native lazy loading so dask is an optional dependency.
        try:
            import dask  # noqa: F401
            chunks: dict | str = {}
        except ImportError:
            chunks = None
        # Choose an explicit engine so files whose extension xarray does not
        # recognise still open.  GeoClaw dtopo/topo files are often named with
        # domain-specific suffixes (e.g. '.dtt3'), and xarray's extension-based
        # engine guessing raises "did not find a match in any of xarray's ...
        # IO backends" for those even though the bytes are valid NetCDF.
        # Prefer netcdf4, fall back to h5netcdf, else leave as None so xarray
        # guesses (and raises its own clear error if no engine is installed).
        engine = next(
            (name for name in ("netcdf4", "h5netcdf")
             if importlib.util.find_spec(name) is not None),
            None,
        )
        # mask_and_scale=True (default) so fill values appear as NaN.
        # decode_timedelta=False: a numeric time coordinate with a bare
        # duration-style units attribute (e.g. "seconds", as written by
        # DTopography._write_netcdf) must stay a plain float array here.
        # Some xarray versions otherwise auto-decode it into timedelta64[ns],
        # and naively casting that to float (as a "seconds" value) reads back
        # the raw nanosecond tick count -- a billion-fold error.
        self.ds: xr.Dataset = xr.open_dataset(
            self.path, engine=engine, chunks=chunks, mask_and_scale=True,
            decode_timedelta=False,
        )

    # ------------------------------------------------------------------
    # Coordinate discovery
    # ------------------------------------------------------------------

    def _find_coord_name(
        self,
        axis_letter: str,
        std_names: list[str],
        fallback_names: list[str],
    ) -> Optional[str]:
        """Find a coordinate by CF conventions: axis > standard_name > name."""
        # 1. Check axis attribute (highest CF priority)
        for name, coord in self.ds.coords.items():
            if coord.attrs.get('axis', '').upper() == axis_letter:
                return str(name)
        # 2. Check standard_name
        for name, coord in self.ds.coords.items():
            if coord.attrs.get('standard_name', '') in std_names:
                return str(name)
        # 3. Fallback: common names (case-insensitive)
        coord_lower = {str(n).lower(): str(n) for n in self.ds.coords}
        for fname in fallback_names:
            if fname.lower() in coord_lower:
                return coord_lower[fname.lower()]
        return None

    def _require_coord(
        self,
        axis_letter: str,
        std_names: list[str],
        fallback_names: list[str],
        label: str,
    ) -> str:
        name = self._find_coord_name(axis_letter, std_names, fallback_names)
        if name is None:
            raise ValueError(
                f"Cannot find {label} coordinate in '{self.path}'.  "
                f"Expected axis='{axis_letter}', standard_name in {std_names}, "
                f"or a name in {fallback_names}.  "
                f"Available coordinates: {list(self.ds.coords)}"
            )
        return name

    def _find_x_name(self) -> str:
        return self._require_coord(
            'X',
            ['longitude'],
            ['longitude', 'lon', 'x', 'nav_lon', 'LONGITUDE', 'LON'],
            'longitude',
        )

    def _find_y_name(self) -> str:
        return self._require_coord(
            'Y',
            ['latitude'],
            ['latitude', 'lat', 'y', 'nav_lat', 'LATITUDE', 'LAT'],
            'latitude',
        )

    def _find_time_name(self) -> Optional[str]:
        return self._find_coord_name(
            'T',
            ['time'],
            ['time', 't', 'TIME', 'valid_time'],
        )

    # ------------------------------------------------------------------
    # Convention detection
    # ------------------------------------------------------------------

    def _detect_lon_wrap(self, x_name: str) -> Optional[int]:
        """Longitude wrap convention for a geographic x axis, else None.

        Returns 360 if longitudes span [0, 360], 180 if [-180, 180].

        The axis is assumed geographic unless there is positive evidence it is
        a projected / rectilinear grid (length units such as ``m``/``km``, a
        projection ``standard_name``, or coordinate values outside the
        plausible degree range [-360, 360]).  In that case it returns None and
        the 0–360 wrap machinery is skipped entirely — such an axis never
        wraps.  Defaulting to geographic preserves behavior for the many files
        whose lon/lat coordinates carry no explicit CF units.
        """
        coord = self.ds[x_name]
        units = str(coord.attrs.get('units', '')).strip().lower()
        std_name = coord.attrs.get('standard_name', '')
        x_min = float(coord.min())
        x_max = float(coord.max())
        projected = (
            units in _PROJECTED_LENGTH_UNITS
            or std_name in _PROJECTED_STANDARD_NAMES
            or x_min < -360.0 - 1e-6
            or x_max > 360.0 + 1e-6
        )
        if projected:
            return None
        return 360 if x_max > 180.0 else 180

    def _detect_y_increasing(self, y_name: str) -> bool:
        """Return True if y increases with array index, else False."""
        y_vals = self.ds[y_name].values
        return bool(y_vals[0] < y_vals[-1])

    def _detect_dim_order(
        self,
        var_name: str,
        x_name: str,
        y_name: str,
        time_name: Optional[str],
    ) -> list[str]:
        """
        Return dimension order for *var_name* using canonical role names
        ('x', 'y', 'time').  Unknown dims are passed through as-is.
        """
        dim_map: dict[str, str] = {x_name: 'x', y_name: 'y'}
        if time_name is not None:
            dim_map[time_name] = 'time'
        return [dim_map.get(d, d) for d in self.ds[var_name].dims]

    # ------------------------------------------------------------------
    # Fill value
    # ------------------------------------------------------------------

    def _resolve_fill_value(self, var_name: str) -> Optional[float]:
        """
        Return fill value with CF-correct precedence: _FillValue > missing_value.

        Note: xarray with mask_and_scale=True has already decoded fill values
        to NaN in the data arrays, but the original attribute is preserved.
        """
        attrs = self.ds[var_name].attrs
        if '_FillValue' in attrs:
            return float(attrs['_FillValue'])
        if 'missing_value' in attrs:
            return float(attrs['missing_value'])
        # xarray may have moved _FillValue to encoding dict
        enc = self.ds[var_name].encoding
        if '_FillValue' in enc:
            return float(enc['_FillValue'])
        return None

    # ------------------------------------------------------------------
    # Crop bounds validation
    # ------------------------------------------------------------------

    def _validate_crop_bounds(
        self,
        x_name: str,
        y_name: str,
        crop_bounds: tuple[float, float, float, float],
    ) -> None:
        """Raise ValueError if crop_bounds exceed file spatial extent."""
        lon_min, lon_max, lat_min, lat_max = crop_bounds
        fx_min = float(self.ds[x_name].min())
        fx_max = float(self.ds[x_name].max())
        fy_min = float(self.ds[y_name].min())
        fy_max = float(self.ds[y_name].max())

        if lon_min < fx_min or lon_max > fx_max:
            raise ValueError(
                f"crop_bounds lon [{lon_min}, {lon_max}] exceed file extent "
                f"[{fx_min}, {fx_max}] in '{self.path}'."
            )
        if lat_min < fy_min or lat_max > fy_max:
            raise ValueError(
                f"crop_bounds lat [{lat_min}, {lat_max}] exceed file extent "
                f"[{fy_min}, {fy_max}] in '{self.path}'."
            )

    # ------------------------------------------------------------------
    # Core inspection
    # ------------------------------------------------------------------

    def inspect(
        self,
        var_name: str,
        time_name_override: Optional[str] = None,
    ) -> FileMetadata:
        """
        Interrogate *var_name* and return a FileMetadata instance.

        Parameters
        ----------
        var_name : str
            Name of the primary data variable (used to determine dim_order).
        time_name_override : str, optional
            Force a specific time coordinate name instead of auto-detecting.
        """
        if var_name not in self.ds:
            raise KeyError(
                f"Variable '{var_name}' not found in '{self.path}'.  "
                f"Available: {list(self.ds.data_vars)}"
            )

        x_name = self._find_x_name()
        y_name = self._find_y_name()
        time_name = (
            time_name_override
            if time_name_override is not None
            else self._find_time_name()
        )

        lon_wrap = self._detect_lon_wrap(x_name)
        y_increasing = self._detect_y_increasing(y_name)
        dim_order = self._detect_dim_order(var_name, x_name, y_name, time_name)
        fill_value = self._resolve_fill_value(var_name)

        if self.crop_bounds is not None:
            self._validate_crop_bounds(x_name, y_name, self.crop_bounds)

        return FileMetadata(
            source_file=self.path,
            x_name=x_name,
            y_name=y_name,
            time_name=time_name,
            lon_wrap=lon_wrap,
            y_increasing=y_increasing,
            dim_order=dim_order,
            fill_value=fill_value,
            crop_bounds=self.crop_bounds,
        )

    def close(self) -> None:
        self.ds.close()

    def __enter__(self) -> 'NetCDFInspector':
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


# ---------------------------------------------------------------------------
# Topo inspector
# ---------------------------------------------------------------------------

class TopoInspector(NetCDFInspector):
    """
    Interrogate a NetCDF bathymetry/topography file.

    Additional checks beyond the base class:

    * Verifies the data variable's units attribute matches the contract unit
      (meters).  If units are convertible via units.py, records the
      source units in the metadata; Fortran will need a conversion factor.
      If units are unrecognised, raises ValueError.
    * Checks for fill values (NaN) within the crop region and raises
      ValueError — silent NaN in bathymetry is numerically fatal.

    Parameters
    ----------
    path : str or Path
        Path to the NetCDF topography file.
    var_name : str, optional
        Name of the elevation/bathymetry variable (e.g. ``'z'``,
        ``'elevation'``).  If omitted, the variable is auto-detected by
        searching for known CF ``standard_name`` values
        (``surface_altitude``, ``height_above_mean_sea_level``, …) and
        then falling back to common names (``z``, ``elevation``, ``topo``,
        …).  A ``ValueError`` is raised if no match is found.
    crop_bounds : tuple of four floats, optional
        (lon_min, lon_max, lat_min, lat_max).
    assume_units : str, optional
        Explicit escape hatch for a file whose elevation variable has *no*
        ``units`` attribute.  When given (e.g. ``'m'``), that unit is assumed
        instead of raising.  This must be set deliberately by the caller;
        missing units are never silently assumed.
    """

    def __init__(
        self,
        path: str | Path,
        var_name: Optional[str] = None,
        crop_bounds: Optional[tuple[float, float, float, float]] = None,
        assume_units: Optional[str] = None,
    ) -> None:
        super().__init__(path, crop_bounds)
        self.var_name = var_name
        self.assume_units = assume_units

    # ------------------------------------------------------------------
    # Variable auto-detection
    # ------------------------------------------------------------------

    #: CF standard_names that indicate an elevation / bathymetry variable,
    #: ordered from most specific to most generic.
    _TOPO_CF_STANDARD_NAMES: tuple[str, ...] = (
        'surface_altitude',
        'height_above_mean_sea_level',
        'height_above_reference_ellipsoid',
        'bedrock_altitude',
        'altitude',
        'height',
        'sea_floor_depth_below_geoid',
    )

    #: Common variable names for elevation / bathymetry data.
    _TOPO_FALLBACK_NAMES: tuple[str, ...] = (
        'z', 'Z',
        'elevation', 'Elevation', 'ELEVATION',
        'topo', 'TOPO',
        'height', 'HEIGHT',
        'altitude', 'ALTITUDE',
        'depth', 'DEPTH',
        'Band1',
        'dem', 'DEM',
        'topography', 'bathymetry', 'bathy',
        'bedrock',
    )

    def _find_topo_var_name(self) -> str:
        """
        Auto-detect the elevation/bathymetry variable.

        Search order:
        1. CF ``standard_name`` attribute against known elevation standard names.
        2. Variable name against ``_TOPO_FALLBACK_NAMES`` (case-sensitive).

        Returns the matched variable name.  Raises ValueError if nothing is found.
        """
        # 1. CF standard_name lookup across all data variables
        for std_name in self._TOPO_CF_STANDARD_NAMES:
            for var_name, var in self.ds.data_vars.items():
                if var.attrs.get('standard_name', '') == std_name:
                    return str(var_name)

        # 2. Name-based fallback
        available = {str(k) for k in self.ds.data_vars}
        for name in self._TOPO_FALLBACK_NAMES:
            if name in available:
                return name

        raise ValueError(
            f"Cannot auto-detect elevation variable in '{self.path}'.  "
            f"No data variable matched known CF standard_names "
            f"{self._TOPO_CF_STANDARD_NAMES} or fallback names "
            f"{self._TOPO_FALLBACK_NAMES}.  "
            f"Pass var_name explicitly.  "
            f"Available data variables: {sorted(available)}"
        )

    # ------------------------------------------------------------------
    # Unit verification
    # ------------------------------------------------------------------

    def _check_topo_units(self) -> str:
        """
        Verify units of *var_name* are meters (GEOCLAW_NETCDF_UNITS['topo']).

        Returns the units string found in the file.  Units are never silently
        assumed or mixed: a missing ``units`` attribute raises ValueError
        (unless *assume_units* was set explicitly), and any non-meter unit --
        even a recognised, convertible one such as ``km`` -- also raises,
        because the Fortran read path (and the in-memory read here) consume the
        raw values as meters with no conversion.  Pre-convert such files to
        meters first.
        """
        contract = GEOCLAW_NETCDF_UNITS['topo']  # 'm'
        cf_unit = self.ds[self.var_name].attrs.get('units', '')

        if not cf_unit:
            if self.assume_units is not None:
                # Deliberate caller override for a known-unitless file.
                return self.assume_units
            raise ValueError(
                f"Variable '{self.var_name}' in '{self.path}' has no 'units' "
                f"attribute.  Units are required and never assumed: add a CF "
                f"'units' attribute (e.g. 'm') to the file, or pass "
                f"assume_units='m' if you are certain the data are in meters."
            )

        if _unit_matches_contract(cf_unit, contract):
            return cf_unit

        # Any other unit -- recognised-but-non-meter (e.g. 'km') or entirely
        # unrecognised -- is rejected: the read paths do not convert, so
        # proceeding would silently corrupt the bathymetry.
        canonical = _normalize_cf_unit(cf_unit)
        if canonical is not None:
            raise ValueError(
                f"Variable '{self.var_name}' in '{self.path}' has units "
                f"'{cf_unit}'; contract requires '{contract}'.  Unit conversion "
                f"on the read path is not supported, so this file is rejected "
                f"rather than silently misread.  Pre-convert the data to meters "
                f"first."
            )
        raise ValueError(
            f"Unrecognised units '{cf_unit}' on variable '{self.var_name}' "
            f"in '{self.path}'.  Contract requires '{contract}'.  Pre-convert "
            f"the file to meters."
        )

    # ------------------------------------------------------------------
    # Fill value check within crop region
    # ------------------------------------------------------------------

    def _check_fill_in_crop(
        self,
        x_name: str,
        y_name: str,
        y_increasing: bool,
    ) -> None:
        """
        Check for NaN / fill values within the crop region.

        This triggers a bounded Dask computation (scalar boolean result).
        Raises ValueError if any fill values are present — silent NaN in
        bathymetry would silently corrupt simulation output.
        """
        var = self.ds[self.var_name]

        if self.crop_bounds is not None:
            lon_min, lon_max, lat_min, lat_max = self.crop_bounds
            # For a decreasing y axis the slice must be reversed.
            lat_slice = (
                slice(lat_min, lat_max)
                if y_increasing
                else slice(lat_max, lat_min)
            )
            var = var.sel({
                x_name: slice(lon_min, lon_max),
                y_name: lat_slice,
            })

        # This .compute() is intentional: we must confirm no fill values exist.
        has_fill = bool(var.isnull().any().compute())
        if has_fill:
            region_str = (
                f"crop region {self.crop_bounds}" if self.crop_bounds is not None
                else "full file extent"
            )
            raise ValueError(
                f"Fill values (NaN) found in topography variable '{self.var_name}' "
                f"within {region_str} of '{self.path}'.  "
                f"Silent NaN in bathymetry would corrupt simulation results.  "
                f"Fill holes before using this file."
            )

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def inspect_topo(self) -> TopoMetadata:
        """
        Fully inspect the topo file and return a TopoMetadata instance.
        """
        if self.var_name is None:
            self.var_name = self._find_topo_var_name()
        base = self.inspect(self.var_name)
        source_units = self._check_topo_units()
        self._check_fill_in_crop(base.x_name, base.y_name, base.y_increasing)

        return TopoMetadata(
            **dataclasses.asdict(base),
            var_name=self.var_name,
            source_units=source_units,
            fill_action='abort',  # fill in topo crop region is always fatal
            lon_wrap_offset=0.0,
        )

    def topo_entries(self) -> list[list]:
        """
        Return a list of ready-to-use topo entries for topofiles.

        Each entry is [4, filepath, TopoMetadata]. When no wrapping is needed
        this returns a list of one entry. When crop_bounds straddle the file's
        lon cut point, returns two entries pointing to the same file with
        different lon_wrap_offset and crop_bounds values.

        crop_bounds on the inspector are in domain coordinates; this method
        converts them to file coordinates before storing in the returned
        metadata. Fortran can then use crop_bounds directly against file
        coordinate arrays before applying lon_wrap_offset.
        """

        # Interrogate without crop validation: self.crop_bounds is in domain
        # coordinates, which may lie outside the file extent.
        saved_crop = self.crop_bounds
        self.crop_bounds = None
        try:
            meta = self.inspect_topo()
        finally:
            self.crop_bounds = saved_crop

        if saved_crop is None:
            return [[4, self.path, dataclasses.replace(meta, lon_wrap_offset=0.0)]]

        assert saved_crop is not None  # narrowing hint: already returned above
        lon_coords = self.ds[meta.x_name].values
        file_lon_min = float(lon_coords.min())
        file_lon_max = float(lon_coords.max())
        if len(lon_coords) > 1:
            lon_resolution = float(abs(lon_coords[1] - lon_coords[0]))
        else:
            lon_resolution = 1e-10  # single-point file, no gap tolerance needed
        crop_lon_min, crop_lon_max, crop_lat_min, crop_lat_max = saved_crop

        # The +/-360 wrap candidates only make sense for a geographic
        # longitude axis.  For a non-geographic x axis (lon_wrap is None,
        # e.g. projected meters) restrict to the identity offset so the crop
        # is taken straight from the file extent.
        entries_spec = _compute_lon_entries(
            file_lon_min, file_lon_max, crop_lon_min, crop_lon_max,
            max_gap=lon_resolution,
            allow_wrap=meta.lon_wrap is not None,
        )

        result = []
        for file_crop_min, file_crop_max, lon_offset in entries_spec:
            new_meta = dataclasses.replace(
                meta,
                crop_bounds=(file_crop_min, file_crop_max, crop_lat_min, crop_lat_max),
                lon_wrap_offset=lon_offset,
            )
            result.append([4, self.path, new_meta])

        return result


# ---------------------------------------------------------------------------
# Met inspector
# ---------------------------------------------------------------------------

class DTopoInspector(NetCDFInspector):
    """
    Interrogate a NetCDF moving-topography (dtopo) file.

    Required roles: longitude, latitude, time, and the deformation (dZ)
    variable.  Beyond the base class this:

    * Discovers the deformation variable (by common names, else the unique
      3-D data variable).
    * Requires time to be the variable's slowest (first) dimension.
    * Validates that the time axis is uniformly spaced — the Fortran dtopo
      machinery reconstructs times as t0 + k*dt.
    * Converts the time axis to simulation seconds: numeric time values pass
      through unchanged; CF datetime values require an explicit
      *time_reference* and become seconds since that reference.

    Parameters
    ----------
    path : str or Path
        Path to the NetCDF file.
    var_name : str, optional
        Name of the deformation variable; auto-discovered if omitted.
    time_reference : str or datetime-like, optional
        Reference datetime when the file's time axis decodes to datetimes.
    crop_bounds : tuple of four floats, optional
        (lon_min, lon_max, lat_min, lat_max); validated against file extent.
    """

    _VAR_NAME_CANDIDATES = ('dz', 'dZ', 'deformation', 'displacement',
                            'dtopo')

    def __init__(
        self,
        path: str | Path,
        var_name: Optional[str] = None,
        time_reference: Optional[str] = None,
        crop_bounds: Optional[tuple[float, float, float, float]] = None,
    ) -> None:
        super().__init__(path, crop_bounds)
        self.var_name = var_name
        self.time_reference = time_reference

    def _find_dtopo_var_name(self) -> str:
        """Find the deformation variable: common names, else unique 3-D var."""
        names_lower = {str(n).lower(): str(n) for n in self.ds.data_vars}
        for cand in self._VAR_NAME_CANDIDATES:
            if cand.lower() in names_lower:
                return names_lower[cand.lower()]

        matches = [str(n) for n, v in self.ds.data_vars.items()
                   if v.ndim == 3]
        if len(matches) == 1:
            return matches[0]
        raise ValueError(
            f"Cannot identify the deformation variable in '{self.path}'.  "
            f"Expected a variable named one of {self._VAR_NAME_CANDIDATES} "
            f"or a unique 3-D data variable; available: "
            f"{list(self.ds.data_vars)}.  Pass var_name explicitly."
        )

    def _compute_time_axis(self, time_name: str) -> tuple[float, float, int]:
        """Return (t0, dt, mt) in simulation seconds; validate uniform dt."""
        tvals = self.ds[time_name].values
        mt = len(tvals)

        if np.issubdtype(tvals.dtype, np.datetime64):
            if self.time_reference is None:
                raise ValueError(
                    f"Time axis '{time_name}' in '{self.path}' decodes to "
                    f"datetimes; pass time_reference to define t=0 in "
                    f"simulation seconds."
                )
            import pandas as pd
            ref = np.datetime64(pd.Timestamp(self.time_reference))
            seconds = (tvals - ref) / np.timedelta64(1, 's')
        elif np.issubdtype(tvals.dtype, np.timedelta64):
            # A duration-style units attribute (e.g. "seconds") decoded to
            # timedelta64 despite decode_timedelta=False (older xarray, or a
            # third-party file with an explicit timedelta64 dtype attribute).
            # Convert properly rather than casting the raw tick count.
            seconds = tvals / np.timedelta64(1, 's')
        else:
            # Numeric time axis (decode_timedelta=False keeps a bare
            # duration-style axis numeric).  Scale to seconds using the
            # coordinate's own 'units' attribute: a "hours"/"minutes" axis
            # would otherwise be read 3600x/60x too fast.  Only when no units
            # attribute is present do we fall back to assuming seconds -- the
            # long-standing contract for a bare numeric dtopo time axis.
            arr = np.asarray(tvals, dtype=float)
            cf_unit = str(self.ds[time_name].attrs.get('units', '')).strip()
            if cf_unit:
                seconds = arr * _cf_time_units_to_seconds_factor(cf_unit)
            else:
                seconds = arr

        if mt < 2:
            return float(seconds[0]), 0.0, mt

        steps = np.diff(seconds)
        dt = float(steps[0])
        if dt <= 0.0 or not np.allclose(steps, dt, rtol=1e-6, atol=0.0):
            raise ValueError(
                f"Time axis '{time_name}' in '{self.path}' is not uniformly "
                f"increasing (steps range [{steps.min()}, {steps.max()}] s). "
                f"The GeoClaw dtopo reader requires uniform time spacing; "
                f"resample the file first."
            )
        return float(seconds[0]), dt, mt

    def inspect_dtopo(self) -> DTopoMetadata:
        """Fully inspect the dtopo file and return a DTopoMetadata."""
        if self.var_name is None:
            self.var_name = self._find_dtopo_var_name()
        base = self.inspect(self.var_name)

        if base.time_name is None:
            raise ValueError(
                f"No time coordinate found in '{self.path}'.  "
                f"dtopography requires a time axis."
            )
        if base.dim_order and base.dim_order[0] != 'time':
            raise ValueError(
                f"Variable '{self.var_name}' in '{self.path}' has dimension "
                f"order {base.dim_order}; the GeoClaw dtopo reader requires "
                f"time as the slowest (first) dimension."
            )

        t0, dt, mt = self._compute_time_axis(base.time_name)

        if base.fill_value is not None:
            warnings.warn(
                f"Deformation variable '{self.var_name}' in '{self.path}' "
                f"declares a fill value ({base.fill_value}); dtopography "
                f"must not contain missing cells.",
                UserWarning,
            )

        return DTopoMetadata(
            **dataclasses.asdict(base),
            var_name=self.var_name,
            t0=t0,
            dt=dt,
            mt=mt,
        )


class MetInspector(NetCDFInspector):
    """
    Interrogate a NetCDF meteorological forcing file.

    Additional checks beyond the base class:

    * Verifies that all requested variables share the same spatial grid and
      time axis.
    * Validates and (if needed) converts units to contract units via units.py.
    * Converts CF time to seconds from a user-provided reference offset.
    * Detects ensemble/member dimensions and raises if any are non-singleton.

    Parameters
    ----------
    path : str or Path
        Path to the NetCDF file.
    variable_map : dict, optional
        Maps GeoClaw role strings to variable names in the file, e.g.::

            {'wind_u': 'u10', 'wind_v': 'v10', 'pressure': 'msl'}

        Roles not given (or the whole map, when omitted) are auto-discovered
        by CF ``standard_name`` first, then by common variable names.  An
        explicit entry overrides discovery for that role.
    crop_bounds : tuple of four floats, optional
        (lon_min, lon_max, lat_min, lat_max).
    time_reference : str or datetime-like, optional
        Reference datetime for the time_offset calculation.  The
        ``time_offset`` written to the descriptor will be the number of
        seconds between *time_reference* and the first time in the file.
        Defaults to the Unix epoch (1970-01-01T00:00:00).
    fill_action : str, optional
        'abort' or 'warn'.  Met files may legitimately have fill values at
        the domain edges (Fortran handles edge fill).  Default is 'warn'.
    assume_units : bool, optional
        Explicit escape hatch for a file whose forcing variables have *no*
        ``units`` attribute.  When True, each variable is assumed to already
        be in its contract unit instead of raising.  This must be set
        deliberately; missing units are never silently assumed.
    """

    # Role discovery, in priority order: CF standard_name, then common
    # variable names (case-insensitive; superset of the historical
    # util.get_netcdf_names fallback lists).
    _MET_STANDARD_NAMES = {
        'wind_u': ['eastward_wind', 'x_wind'],
        'wind_v': ['northward_wind', 'y_wind'],
        'pressure': ['air_pressure_at_mean_sea_level',
                     'air_pressure_at_sea_level',
                     'surface_air_pressure'],
    }
    _MET_FALLBACK_NAMES = {
        'wind_u': ['wind_x', 'wind_u', 'u10', 'u'],
        'wind_v': ['wind_y', 'wind_v', 'v10', 'v'],
        'pressure': ['pressure', 'msl', 'sp'],
    }

    def __init__(
        self,
        path: str | Path,
        variable_map: Optional[dict[str, str]] = None,
        crop_bounds: Optional[tuple[float, float, float, float]] = None,
        time_reference: Optional[str] = None,
        fill_action: str = 'warn',
        assume_units: bool = False,
    ) -> None:
        super().__init__(path, crop_bounds)
        # {geoclaw_role: var_name_in_file}; missing roles auto-discovered
        # in inspect_met() via _resolve_variable_map()
        self.variable_map = variable_map
        self.time_reference = time_reference
        if fill_action not in ('abort', 'warn'):
            raise ValueError(f"fill_action must be 'abort' or 'warn', got '{fill_action}'")
        self.fill_action = fill_action
        self.assume_units = assume_units

    # ------------------------------------------------------------------
    # Variable role discovery
    # ------------------------------------------------------------------

    def _resolve_variable_map(
        self,
        overrides: Optional[dict[str, str]],
    ) -> dict[str, str]:
        """Return a complete {role: var_name} map.

        Explicit *overrides* entries are validated against the file and used
        as-is; the remaining roles are discovered by CF ``standard_name``
        first, then by common variable names.
        """
        overrides = dict(overrides or {})
        available = {str(n) for n in self.ds.data_vars} | \
                    {str(n) for n in self.ds.coords}

        bad = {role: name for role, name in overrides.items()
               if name not in available}
        if bad:
            raise ValueError(
                f"The following names in variable_map were not found in "
                f"'{self.path}':\n"
                + "\n".join(f"  '{role}': '{name}'"
                            for role, name in bad.items())
                + f"\nAvailable variables: {sorted(available)}"
            )

        result = dict(overrides)
        names_lower = {str(n).lower(): str(n) for n in self.ds.data_vars}
        for role in ('wind_u', 'wind_v', 'pressure'):
            if role in result:
                continue
            # 1. CF standard_name
            for name, var in self.ds.data_vars.items():
                if var.attrs.get('standard_name', '') in \
                        self._MET_STANDARD_NAMES[role]:
                    result[role] = str(name)
                    break
            if role in result:
                continue
            # 2. Common names (case-insensitive)
            for cand in self._MET_FALLBACK_NAMES[role]:
                if cand in names_lower:
                    result[role] = names_lower[cand]
                    break
            if role not in result:
                raise ValueError(
                    f"Cannot find a variable for met role '{role}' in "
                    f"'{self.path}'.  Expected standard_name in "
                    f"{self._MET_STANDARD_NAMES[role]} or a name in "
                    f"{self._MET_FALLBACK_NAMES[role]}; available: "
                    f"{list(self.ds.data_vars)}.  Pass variable_map to "
                    f"specify it explicitly."
                )
        return result

    # ------------------------------------------------------------------
    # Variable grid/time consistency
    # ------------------------------------------------------------------

    def _check_variable_consistency(
        self,
        x_name: str,
        y_name: str,
        time_name: str,
    ) -> None:
        """
        Verify all variables share the same spatial grid and time axis.

        Raises ValueError on any mismatch.
        """
        ref_role = next(iter(self.variable_map))
        ref_var_name = self.variable_map[ref_role]
        ref_dims = set(self.ds[ref_var_name].dims)

        for role, var_name in self.variable_map.items():
            if var_name not in self.ds:
                raise KeyError(
                    f"Variable '{var_name}' (role '{role}') not found in "
                    f"'{self.path}'.  Available: {list(self.ds.data_vars)}"
                )
            dims = set(self.ds[var_name].dims)
            if dims != ref_dims:
                raise ValueError(
                    f"Variable '{var_name}' (role '{role}') has dims {dims} "
                    f"but '{ref_var_name}' (role '{ref_role}') has dims {ref_dims}.  "
                    f"All met variables must share the same grid."
                )
            # Verify that the expected coordinate dims are actually present
            for coord_name, label in [
                (x_name, 'longitude'),
                (y_name, 'latitude'),
                (time_name, 'time'),
            ]:
                if coord_name not in dims:
                    raise ValueError(
                        f"Variable '{var_name}' dims {dims} do not include the "
                        f"{label} coordinate '{coord_name}'."
                    )

    # ------------------------------------------------------------------
    # Ensemble dimension check
    # ------------------------------------------------------------------

    def _check_ensemble_dims(
        self,
        known_dims: set[str],
    ) -> None:
        """
        Raise ValueError if any variable has non-singleton extra dimensions.

        Extra dimensions (beyond lon, lat, time) are assumed to be ensemble /
        member axes, which GeoClaw does not support.
        """
        for role, var_name in self.variable_map.items():
            var = self.ds[var_name]
            extra = [d for d in var.dims if d not in known_dims]
            for dim in extra:
                size = self.ds.sizes[dim]
                if size != 1:
                    raise ValueError(
                        f"Variable '{var_name}' (role '{role}') has non-singleton "
                        f"dimension '{dim}' (size {size}).  Ensemble / member "
                        f"dimensions are not supported by GeoClaw.  "
                        f"Select a single member before inspecting."
                    )

    # ------------------------------------------------------------------
    # Unit verification / conversion
    # ------------------------------------------------------------------

    def _check_met_units(self) -> list[MetVariableInfo]:
        """
        Verify units for each variable match its contract unit.

        Returns a list of MetVariableInfo with source_units populated.  Units
        are never silently assumed or mixed: a missing ``units`` attribute
        raises ValueError (unless *assume_units* was set explicitly), and any
        non-contract unit -- even a recognised, convertible one such as
        ``hPa`` or ``knots`` -- also raises, because the Fortran read path
        consumes the raw values with no conversion.  Pre-convert such files to
        the contract units first.
        """
        result: list[MetVariableInfo] = []
        for role, var_name in self.variable_map.items():
            contract = GEOCLAW_NETCDF_UNITS.get(role)
            if contract is None:
                # Role not in contract (e.g. future extension) — just record
                warnings.warn(
                    f"No contract unit defined for role '{role}'.  "
                    f"Skipping unit check for '{var_name}'.",
                    stacklevel=3,
                )
                source_units = self.ds[var_name].attrs.get('units', 'unknown')
                result.append(MetVariableInfo(
                    var_name=var_name,
                    geoclaw_role=role,
                    source_units=source_units,
                ))
                continue

            cf_unit = self.ds[var_name].attrs.get('units', '')
            if not cf_unit:
                if self.assume_units:
                    # Deliberate caller override for a known-unitless file.
                    result.append(MetVariableInfo(
                        var_name=var_name,
                        geoclaw_role=role,
                        source_units=contract,
                    ))
                    continue
                raise ValueError(
                    f"Variable '{var_name}' (role '{role}') in '{self.path}' "
                    f"has no 'units' attribute.  Units are required and never "
                    f"assumed: add a CF 'units' attribute (e.g. '{contract}') "
                    f"to the file, or pass assume_units=True if you are certain "
                    f"the data are already in '{contract}'."
                )

            if _unit_matches_contract(cf_unit, contract):
                result.append(MetVariableInfo(
                    var_name=var_name,
                    geoclaw_role=role,
                    source_units=cf_unit,
                ))
                continue

            # Any other unit -- recognised-but-non-contract (e.g. 'hPa',
            # 'knots') or entirely unrecognised -- is rejected: the read path
            # does not convert, so proceeding would silently corrupt the
            # forcing.  Pre-convert the file to the contract units first.
            canonical = _normalize_cf_unit(cf_unit)
            if canonical is not None:
                raise ValueError(
                    f"Variable '{var_name}' (role '{role}') in '{self.path}' "
                    f"has units '{cf_unit}'; contract requires '{contract}'.  "
                    f"Unit conversion on the read path is not supported, so "
                    f"this file is rejected rather than silently misread.  "
                    f"Pre-convert the data to '{contract}' first."
                )
            raise ValueError(
                f"Unrecognised units '{cf_unit}' on variable '{var_name}' "
                f"(role '{role}') in '{self.path}'.  Contract requires "
                f"'{contract}'.  Pre-convert the file to '{contract}'."
            )

        return result

    # ------------------------------------------------------------------
    # CF time -> seconds offset
    # ------------------------------------------------------------------

    def _compute_time_offset(self, time_name: str) -> float:
        """
        Compute time_offset in seconds.

        time_offset = seconds between *time_reference* and the first time
        value in the file (after CF decoding by xarray).

        This reads only the first element of the time coordinate — a minimal
        data fetch, not a full array load.
        """
        import pandas as pd

        time_coord = self.ds[time_name]
        # Access only the first element (.values[0] fetches one chunk element)
        t0_raw = time_coord.values[0]

        # A bare numeric time axis (e.g. units "hours" with no reference date)
        # is not an absolute time.  Feeding it to pd.Timestamp would silently
        # interpret the raw number as nanoseconds and produce a wildly wrong
        # offset (the same class of bug as the dtopo time-axis reader), so
        # reject it explicitly rather than assuming.  A CF datetime axis
        # ("<unit> since <date>") is decoded by xarray to datetime64 and passes
        # through; cftime objects (non-standard calendars) remain object dtype
        # and are handled by the try/except below.
        _t0_arr = np.asarray(t0_raw)
        if (np.issubdtype(_t0_arr.dtype, np.floating)
                or np.issubdtype(_t0_arr.dtype, np.integer)):
            cf_unit = str(time_coord.attrs.get('units', '')).strip()
            raise ValueError(
                f"Time axis '{time_name}' in '{self.path}' is numeric "
                f"(units {cf_unit!r}); met/storm forcing needs an absolute time "
                f"reference.  Give the time coordinate CF units of the form "
                f"'seconds since <date>' (e.g. 'seconds since 2020-01-01'), not "
                f"a bare duration."
            )

        # The met time axis must be in seconds.  Fortran reads the raw time
        # values into an integer array and computes raw[i] - raw[0] +
        # nc_time_offset, treating them as seconds without scaling by the CF
        # unit, so an "hours since"/"days since" axis (e.g. a raw ERA5 file)
        # would be silently misread.  xarray moves the original units string to
        # .encoding when it decodes the datetime axis; reject any non-second
        # unit here rather than producing wrong forcing.
        _raw_units = str(time_coord.encoding.get('units')
                         or time_coord.attrs.get('units', '')).strip()
        if ' since ' in _raw_units.lower():
            _unit_part = _raw_units.lower().split(' since ', 1)[0].strip()
            if _unit_part not in ('s', 'sec', 'secs', 'second', 'seconds'):
                raise ValueError(
                    f"Time axis '{time_name}' in '{self.path}' has units "
                    f"'{_raw_units}'; met/storm forcing requires the time axis "
                    f"in seconds ('seconds since <date>').  GeoClaw reads the "
                    f"raw time values as integer seconds without unit scaling, "
                    f"so convert the file to 'seconds since ...' first."
                )

        try:
            t0 = pd.Timestamp(t0_raw)
        except Exception as exc:
            raise ValueError(
                f"Cannot convert first time value '{t0_raw}' in '{self.path}' "
                f"to a pandas Timestamp.  Ensure xarray can decode CF time "
                f"(check 'units' and 'calendar' attributes on '{time_name}')."
            ) from exc

        if self.time_reference is None:
            ref = pd.Timestamp('1970-01-01T00:00:00')  # Unix epoch
        else:
            ref = pd.Timestamp(self.time_reference)

        offset_s = (t0 - ref).total_seconds()
        return float(offset_s)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def inspect_met(self) -> MetMetadata:
        """
        Fully inspect the met file and return a MetMetadata instance.

        Roles missing from variable_map are auto-discovered (CF
        standard_name first, then common variable names).
        """
        self.variable_map = self._resolve_variable_map(self.variable_map)

        # Use the first variable for base coordinate detection
        primary_var = next(iter(self.variable_map.values()))
        base = self.inspect(primary_var)

        if base.time_name is None:
            raise ValueError(
                f"No time coordinate found in '{self.path}'.  "
                f"Met forcing requires a time axis."
            )

        known_dims = {base.x_name, base.y_name, base.time_name}
        self._check_variable_consistency(
            base.x_name, base.y_name, base.time_name
        )
        self._check_ensemble_dims(known_dims)
        variable_infos = self._check_met_units()
        time_offset = self._compute_time_offset(base.time_name)

        return MetMetadata(
            **dataclasses.asdict(base),
            variables=variable_infos,
            time_offset=time_offset,
            fill_action=self.fill_action,
        )


# ---------------------------------------------------------------------------
# CF Normalizer
# ---------------------------------------------------------------------------

class CFNormalizer:
    """
    Normalize CF metadata in an xarray Dataset.

    Performs in-memory attribute fixups:
    * Renames coordinate variables to CF standard names where unambiguous.
    * Adds ``standard_name``, ``axis``, and ``units`` attributes to coordinate
      variables if missing.
    * Resolves ``_FillValue`` vs ``missing_value`` conflicts (CF precedence:
      ``_FillValue`` wins; ``missing_value`` is promoted if ``_FillValue``
      absent; conflict triggers a warning).

    Does NOT resample, reproject, or modify data values.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to normalise.  A copy is made; the original is not modified.

    Examples
    --------
    >>> ds_norm = CFNormalizer(ds).normalize()
    """

    # Alias sets for unambiguous renaming
    _LON_ALIASES: frozenset[str] = frozenset(
        {'lon', 'longitude', 'x', 'nav_lon', 'LONGITUDE', 'LON', 'Longitude'}
    )
    _LAT_ALIASES: frozenset[str] = frozenset(
        {'lat', 'latitude', 'y', 'nav_lat', 'LATITUDE', 'LAT', 'Latitude'}
    )
    _TIME_ALIASES: frozenset[str] = frozenset(
        {'time', 't', 'TIME', 'valid_time', 'Time'}
    )

    def __init__(self, ds: xr.Dataset) -> None:
        self.ds = ds.copy()

    def normalize(self) -> xr.Dataset:
        """Apply all normalizations and return the modified dataset."""
        self._rename_coords()
        self._fix_coord_attrs()
        self._resolve_fill_values()
        return self.ds

    def _rename_coords(self) -> None:
        """Rename coordinate variables to CF standard names."""
        rename_map: dict[str, str] = {}
        for name in list(self.ds.coords):
            sname = str(name)
            if sname in self._LON_ALIASES and sname != 'longitude':
                if 'longitude' not in self.ds.coords:
                    rename_map[sname] = 'longitude'
            elif sname in self._LAT_ALIASES and sname != 'latitude':
                if 'latitude' not in self.ds.coords:
                    rename_map[sname] = 'latitude'
            elif sname in self._TIME_ALIASES and sname != 'time':
                if 'time' not in self.ds.coords:
                    rename_map[sname] = 'time'
        if rename_map:
            self.ds = self.ds.rename(rename_map)

    def _fix_coord_attrs(self) -> None:
        """Add standard_name, axis, and units to recognized coordinates."""
        _COORD_STANDARDS: dict[str, tuple[str, Optional[str], str]] = {
            'longitude': ('longitude', 'degrees_east',  'X'),
            'latitude':  ('latitude',  'degrees_north', 'Y'),
            'time':      ('time',      None,            'T'),
        }
        for coord_name, (std_name, units_val, axis) in _COORD_STANDARDS.items():
            if coord_name not in self.ds.coords:
                continue
            attrs = dict(self.ds[coord_name].attrs)
            attrs.setdefault('standard_name', std_name)
            attrs.setdefault('axis', axis)
            if units_val is not None:
                attrs.setdefault('units', units_val)
            self.ds[coord_name].attrs = attrs

    def _resolve_fill_values(self) -> None:
        """
        Resolve _FillValue / missing_value conflicts for all data variables.

        CF precedence: _FillValue wins.  If only missing_value is present it
        is promoted to _FillValue.  Conflicting values trigger a UserWarning.
        """
        for var_name in list(self.ds.data_vars):
            attrs = dict(self.ds[var_name].attrs)
            fv: Optional[float] = attrs.get('_FillValue')
            mv: Optional[float] = attrs.get('missing_value')

            if fv is not None and mv is not None:
                if fv != mv:
                    warnings.warn(
                        f"Variable '{var_name}' has conflicting "
                        f"_FillValue={fv} and missing_value={mv}.  "
                        f"Using _FillValue (CF precedence); "
                        f"missing_value removed.",
                        UserWarning,
                        stacklevel=2,
                    )
                del attrs['missing_value']
                self.ds[var_name].attrs = attrs

            elif mv is not None:
                # Promote missing_value -> _FillValue
                attrs['_FillValue'] = mv
                del attrs['missing_value']
                self.ds[var_name].attrs = attrs


# ---------------------------------------------------------------------------
# Longitude entry computation
# ---------------------------------------------------------------------------

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

    lon_offset is the scalar Fortran adds to file coordinates to produce
    domain coordinates: x_domain = x_file + lon_offset.

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


# ---------------------------------------------------------------------------
# Descriptor writer
# ---------------------------------------------------------------------------

class DescriptorWriter:
    """
    Write NetCDF descriptor metadata for GeoClaw input files.

    For topo (type 4) entries the descriptor is a block of ``key = value``
    lines written immediately after the ``topo_type`` line in *topo.data*.
    Fortran's ``read_netcdf_descriptor`` parses these lines until it
    encounters a blank line.

    Usage — topo::

        meta = TopoInspector(path, var_name='z').inspect_topo()
        with open('topo.data', 'a') as f:
            DescriptorWriter.write_topo_descriptor(f, meta)

    Usage — met (storm file)::

        meta = MetInspector(path, var_map).inspect_met()
        with open('storm.nc', 'w') as f:
            f.write('netcdf\\n')   # format header written by caller
            DescriptorWriter.write_met_descriptor(f, meta)
    """

    @staticmethod
    def write_topo_descriptor(f, meta: 'TopoMetadata') -> None:
        """
        Write the key=value descriptor block for one topo file.

        *f* should be positioned immediately after the ``topo_type`` line
        has been written.  A trailing blank line is written to terminate
        the block for the Fortran parser.

        Parameters
        ----------
        f : file-like object
            Open for writing (text mode).
        meta : TopoMetadata
            Output of ``TopoInspector.inspect_topo()``.
        """
        # Fortran adds lon_wrap_offset to all file coordinate values after
        # reading:
        #   x_domain = x_file + lon_wrap_offset
        # crop_bounds written here are in FILE coordinates (converted from
        # domain coordinates by topo_entries()).
        f.write(f"var_name       = {meta.var_name}\n")
        f.write(f"x_name         = {meta.x_name}\n")
        f.write(f"y_name         = {meta.y_name}\n")
        f.write(f"lon_wrap_offset = {meta.lon_wrap_offset!r}\n")
        f.write(f"y_increasing   = {meta.y_increasing}\n")
        f.write(f"dim_order      = {','.join(meta.dim_order)}\n")
        if meta.fill_value is not None:
            f.write(f"fill_value     = {meta.fill_value!r}\n")
        f.write(f"fill_action    = {meta.fill_action}\n")
        if meta.crop_bounds is not None:
            x0, x1, y0, y1 = meta.crop_bounds
            f.write(f"crop_bounds    = {x0} {x1} {y0} {y1}\n")
        f.write("\n")  # blank line terminates block for Fortran parser

    @staticmethod
    def write_dtopo_descriptor(f, meta: 'DTopoMetadata') -> None:
        """
        Write the key=value descriptor block for one dtopo file.

        Same format and parser contract as the topo descriptor (a blank
        line terminates the block), written immediately after the 9-line
        per-file block in *dtopo.data*.  The time axis is carried as
        (t0, dt) in simulation seconds so Fortran never parses CF time.

        Parameters
        ----------
        f : file-like object
            Open for writing (text mode).
        meta : DTopoMetadata
            Output of ``DTopoInspector.inspect_dtopo()``.
        """
        f.write(f"var_name       = {meta.var_name}\n")
        f.write(f"x_name         = {meta.x_name}\n")
        f.write(f"y_name         = {meta.y_name}\n")
        f.write(f"time_name      = {meta.time_name}\n")
        f.write(f"lon_wrap_offset = {meta.lon_wrap_offset!r}\n")
        f.write(f"y_increasing   = {meta.y_increasing}\n")
        f.write(f"dim_order      = {','.join(meta.dim_order)}\n")
        f.write(f"t0             = {meta.t0!r}\n")
        f.write(f"dt             = {meta.dt!r}\n")
        f.write("\n")  # blank line terminates block for Fortran parser

    @staticmethod
    def write_met_descriptor(f, meta: 'MetMetadata') -> None:
        """
        Write the ``&file_info`` / ``&variable_info`` namelist-style body
        of a GeoClaw NetCDF met/storm descriptor file.

        The caller is responsible for writing the ``netcdf`` format header
        line before calling this method.

        Parameters
        ----------
        f : file-like object
            Open for writing (text mode).
        meta : MetMetadata
            Output of ``MetInspector.inspect_met()``.
        """
        f.write("&file_info\n")
        f.write(f"  x_name         = {meta.x_name}\n")
        f.write(f"  y_name         = {meta.y_name}\n")
        f.write(f"  time_name      = {meta.time_name}\n")
        f.write(f"  dim_order      = {','.join(meta.dim_order)}\n")
        # lon_wrap only applies to a geographic longitude axis.  Omit it for a
        # non-geographic x axis (lon_wrap is None): Fortran defaults to 180
        # (no 0-360 normalization), which is the correct no-wrap behavior.
        if meta.lon_wrap is not None:
            f.write(f"  lon_wrap       = {meta.lon_wrap}\n")
        f.write(f"  y_increasing   = {meta.y_increasing}\n")
        if meta.fill_value is not None:
            f.write(f"  fill_value     = {meta.fill_value!r}\n")
        f.write(f"  fill_action    = {meta.fill_action}\n")
        f.write(f"  time_offset    = {meta.time_offset!r}\n")
        if meta.crop_bounds is not None:
            x0, x1, y0, y1 = meta.crop_bounds
            f.write(f"  crop_bounds    = {x0} {x1} {y0} {y1}\n")
        f.write("/\n")
        for var in meta.variables:
            f.write(
                f"&variable_info"
                f"  var_name={var.var_name}"
                f"  geoclaw_role={var.geoclaw_role}"
                f"  /\n"
            )
