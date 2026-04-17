#!/usr/bin/env python
# encoding: utf-8
r"""
NetCDF utilities for GeoClaw NetCDF input.

Provides interrogator classes that inspect NetCDF files for coordinate
metadata, unit consistency, and fill values without loading data arrays.
The interrogators produce metadata dataclasses consumed by DescriptorWriter
(not yet implemented — pending confirmation of Fortran topo.data parser).

Classes
-------
NetCDFInterrogator
    Base class: opens a file, discovers coordinate variables, detects lon/lat
    conventions, resolves fill value, validates crop bounds.

TopoInterrogator(NetCDFInterrogator)
    Adds bathymetry-specific checks: fill values within crop region (fatal),
    unit verification against GEOCLAW_NETCDF_UNITS['topo'].

MetInterrogator(NetCDFInterrogator)
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

import dataclasses
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


def _normalize_cf_unit(cf_unit: str) -> Optional[str]:
    """Return the units.py abbreviation for *cf_unit*, or None if unknown."""
    return _CF_TO_UNITS_PY.get(cf_unit.strip())


def _unit_matches_contract(cf_unit: str, contract_unit: str) -> bool:
    """Return True if *cf_unit* is equivalent to *contract_unit*."""
    aliases = _CONTRACT_UNIT_CF_ALIASES.get(contract_unit, frozenset())
    return cf_unit.strip() in aliases


# ---------------------------------------------------------------------------
# Metadata dataclasses
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class FileMetadata:
    """Coordinate and convention metadata for a single NetCDF file."""
    source_file: Path
    lon_name: str
    lat_name: str
    time_name: Optional[str]
    lon_convention: int          # 180 -> [-180, 180];  360 -> [0, 360]
    lat_order: str               # 'N_to_S' or 'S_to_N'
    dim_order: list[str]         # canonical role names, e.g. ['lat', 'lon']
    fill_value: Optional[float]  # resolved from _FillValue / missing_value
    crop_bounds: Optional[tuple[float, float, float, float]]  # lon0,lon1,lat0,lat1


@dataclasses.dataclass
class TopoMetadata(FileMetadata):
    """FileMetadata plus topo-specific fields."""
    var_name: str
    source_units: str    # units as found in file
    fill_action: str     # 'abort' (only value for topo; fill = fatal)


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
# Base interrogator
# ---------------------------------------------------------------------------

class NetCDFInterrogator:
    """
    Open a NetCDF file and interrogate its coordinate metadata.

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
        # mask_and_scale=True (default) so fill values appear as NaN.
        self.ds: xr.Dataset = xr.open_dataset(
            self.path, chunks=chunks, mask_and_scale=True
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

    def _find_lon_name(self) -> str:
        return self._require_coord(
            'X',
            ['longitude'],
            ['longitude', 'lon', 'x', 'nav_lon', 'LONGITUDE', 'LON'],
            'longitude',
        )

    def _find_lat_name(self) -> str:
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

    def _detect_lon_convention(self, lon_name: str) -> int:
        """Return 180 if longitudes are in [-180, 180], else 360."""
        lon_max = float(self.ds[lon_name].max())
        return 360 if lon_max > 180.0 else 180

    def _detect_lat_order(self, lat_name: str) -> str:
        """Return 'N_to_S' if latitudes decrease, else 'S_to_N'."""
        lat_vals = self.ds[lat_name].values
        if lat_vals[0] > lat_vals[-1]:
            return 'N_to_S'
        return 'S_to_N'

    def _detect_dim_order(
        self,
        var_name: str,
        lon_name: str,
        lat_name: str,
        time_name: Optional[str],
    ) -> list[str]:
        """
        Return dimension order for *var_name* using canonical role names
        ('lon', 'lat', 'time').  Unknown dims are passed through as-is.
        """
        dim_map: dict[str, str] = {lon_name: 'lon', lat_name: 'lat'}
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
        lon_name: str,
        lat_name: str,
        crop_bounds: tuple[float, float, float, float],
    ) -> None:
        """Raise ValueError if crop_bounds exceed file spatial extent."""
        lon_min, lon_max, lat_min, lat_max = crop_bounds
        flon_min = float(self.ds[lon_name].min())
        flon_max = float(self.ds[lon_name].max())
        flat_min = float(self.ds[lat_name].min())
        flat_max = float(self.ds[lat_name].max())

        if lon_min < flon_min or lon_max > flon_max:
            raise ValueError(
                f"crop_bounds lon [{lon_min}, {lon_max}] exceed file extent "
                f"[{flon_min}, {flon_max}] in '{self.path}'."
            )
        if lat_min < flat_min or lat_max > flat_max:
            raise ValueError(
                f"crop_bounds lat [{lat_min}, {lat_max}] exceed file extent "
                f"[{flat_min}, {flat_max}] in '{self.path}'."
            )

    # ------------------------------------------------------------------
    # Core interrogation
    # ------------------------------------------------------------------

    def interrogate(
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

        lon_name = self._find_lon_name()
        lat_name = self._find_lat_name()
        time_name = (
            time_name_override
            if time_name_override is not None
            else self._find_time_name()
        )

        lon_convention = self._detect_lon_convention(lon_name)
        lat_order = self._detect_lat_order(lat_name)
        dim_order = self._detect_dim_order(var_name, lon_name, lat_name, time_name)
        fill_value = self._resolve_fill_value(var_name)

        if self.crop_bounds is not None:
            self._validate_crop_bounds(lon_name, lat_name, self.crop_bounds)

        return FileMetadata(
            source_file=self.path,
            lon_name=lon_name,
            lat_name=lat_name,
            time_name=time_name,
            lon_convention=lon_convention,
            lat_order=lat_order,
            dim_order=dim_order,
            fill_value=fill_value,
            crop_bounds=self.crop_bounds,
        )

    def close(self) -> None:
        self.ds.close()

    def __enter__(self) -> 'NetCDFInterrogator':
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


# ---------------------------------------------------------------------------
# Topo interrogator
# ---------------------------------------------------------------------------

class TopoInterrogator(NetCDFInterrogator):
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
    var_name : str
        Name of the elevation/bathymetry variable (e.g. 'z', 'elevation').
    crop_bounds : tuple of four floats, optional
        (lon_min, lon_max, lat_min, lat_max).
    """

    def __init__(
        self,
        path: str | Path,
        var_name: str,
        crop_bounds: Optional[tuple[float, float, float, float]] = None,
    ) -> None:
        super().__init__(path, crop_bounds)
        self.var_name = var_name

    # ------------------------------------------------------------------
    # Unit verification
    # ------------------------------------------------------------------

    def _check_topo_units(self) -> str:
        """
        Verify units of *var_name* match GEOCLAW_NETCDF_UNITS['topo'] ('m').

        Returns the units string found in the file.  Raises ValueError if the
        units are not recognisable or not convertible to meters.
        """
        contract = GEOCLAW_NETCDF_UNITS['topo']  # 'm'
        cf_unit = self.ds[self.var_name].attrs.get('units', '')

        if not cf_unit:
            warnings.warn(
                f"Variable '{self.var_name}' in '{self.path}' has no 'units' "
                f"attribute.  Assuming '{contract}' (contract unit for topo).  "
                f"If the data are not in meters, results will be incorrect.",
                stacklevel=3,
            )
            return contract

        if _unit_matches_contract(cf_unit, contract):
            return cf_unit

        # Try to find a conversion path via units.py
        canonical = _normalize_cf_unit(cf_unit)
        if canonical is None:
            raise ValueError(
                f"Unrecognised units '{cf_unit}' on variable '{self.var_name}' "
                f"in '{self.path}'.  Contract requires '{contract}'.  "
                f"Add the unit string to _CF_TO_UNITS_PY or pre-convert the file."
            )
        # Verify units.convert() can handle the conversion (raises if not)
        try:
            units_convert(1.0, canonical, contract)
        except (ValueError, KeyError) as exc:
            raise ValueError(
                f"Cannot convert units '{cf_unit}' -> '{contract}' for variable "
                f"'{self.var_name}' in '{self.path}': {exc}"
            ) from exc

        warnings.warn(
            f"Variable '{self.var_name}' has units '{cf_unit}'; contract is "
            f"'{contract}'.  A unit conversion will be needed before Fortran "
            f"reads this file.  Pre-convert the data or implement unit scaling "
            f"in the descriptor.",
            stacklevel=3,
        )
        return cf_unit

    # ------------------------------------------------------------------
    # Fill value check within crop region
    # ------------------------------------------------------------------

    def _check_fill_in_crop(
        self,
        lon_name: str,
        lat_name: str,
        lat_order: str,
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
            # For a decreasing lat axis (N_to_S) the slice must be reversed.
            lat_slice = (
                slice(lat_max, lat_min)
                if lat_order == 'N_to_S'
                else slice(lat_min, lat_max)
            )
            var = var.sel({
                lon_name: slice(lon_min, lon_max),
                lat_name: lat_slice,
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

    def interrogate_topo(self) -> TopoMetadata:
        """
        Fully interrogate the topo file and return a TopoMetadata instance.
        """
        base = self.interrogate(self.var_name)
        source_units = self._check_topo_units()
        self._check_fill_in_crop(base.lon_name, base.lat_name, base.lat_order)

        return TopoMetadata(
            **dataclasses.asdict(base),
            var_name=self.var_name,
            source_units=source_units,
            fill_action='abort',  # fill in topo crop region is always fatal
        )


# ---------------------------------------------------------------------------
# Met interrogator
# ---------------------------------------------------------------------------

class MetInterrogator(NetCDFInterrogator):
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
    variable_map : dict
        Maps GeoClaw role strings to variable names in the file, e.g.::

            {'wind_u': 'u10', 'wind_v': 'v10', 'pressure': 'msl'}

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
    """

    def __init__(
        self,
        path: str | Path,
        variable_map: dict[str, str],
        crop_bounds: Optional[tuple[float, float, float, float]] = None,
        time_reference: Optional[str] = None,
        fill_action: str = 'warn',
    ) -> None:
        super().__init__(path, crop_bounds)
        self.variable_map = variable_map  # {geoclaw_role: var_name_in_file}
        self.time_reference = time_reference
        if fill_action not in ('abort', 'warn'):
            raise ValueError(f"fill_action must be 'abort' or 'warn', got '{fill_action}'")
        self.fill_action = fill_action

    # ------------------------------------------------------------------
    # Variable grid/time consistency
    # ------------------------------------------------------------------

    def _check_variable_consistency(
        self,
        lon_name: str,
        lat_name: str,
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
                (lon_name, 'longitude'),
                (lat_name, 'latitude'),
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
                        f"Select a single member before interrogating."
                    )

    # ------------------------------------------------------------------
    # Unit verification / conversion
    # ------------------------------------------------------------------

    def _check_met_units(self) -> list[MetVariableInfo]:
        """
        Verify or convert units for each variable to contract units.

        Returns a list of MetVariableInfo with source_units populated.
        Raises ValueError for unrecognised or unconvertible units.
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
                warnings.warn(
                    f"Variable '{var_name}' (role '{role}') has no 'units' "
                    f"attribute.  Assuming contract unit '{contract}'.  "
                    f"Verify the data are actually in '{contract}'.",
                    stacklevel=3,
                )
                result.append(MetVariableInfo(
                    var_name=var_name,
                    geoclaw_role=role,
                    source_units=contract,
                ))
                continue

            if _unit_matches_contract(cf_unit, contract):
                result.append(MetVariableInfo(
                    var_name=var_name,
                    geoclaw_role=role,
                    source_units=cf_unit,
                ))
                continue

            # Attempt conversion path
            canonical = _normalize_cf_unit(cf_unit)
            if canonical is None:
                raise ValueError(
                    f"Unrecognised units '{cf_unit}' on variable '{var_name}' "
                    f"(role '{role}') in '{self.path}'.  "
                    f"Contract requires '{contract}'.  "
                    f"Add the unit to _CF_TO_UNITS_PY or pre-convert the file."
                )
            try:
                units_convert(1.0, canonical, contract)
            except (ValueError, KeyError) as exc:
                raise ValueError(
                    f"Cannot convert units '{cf_unit}' -> '{contract}' for "
                    f"variable '{var_name}' (role '{role}'): {exc}"
                ) from exc

            warnings.warn(
                f"Variable '{var_name}' (role '{role}') has units '{cf_unit}'; "
                f"contract is '{contract}'.  A unit conversion factor must be "
                f"applied before Fortran reads this data.",
                stacklevel=3,
            )
            result.append(MetVariableInfo(
                var_name=var_name,
                geoclaw_role=role,
                source_units=cf_unit,
            ))

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

    def interrogate_met(self) -> MetMetadata:
        """
        Fully interrogate the met file and return a MetMetadata instance.
        """
        if not self.variable_map:
            raise ValueError("variable_map must not be empty.")

        # Use the first variable for base coordinate detection
        primary_var = next(iter(self.variable_map.values()))
        base = self.interrogate(primary_var)

        if base.time_name is None:
            raise ValueError(
                f"No time coordinate found in '{self.path}'.  "
                f"Met forcing requires a time axis."
            )

        known_dims = {base.lon_name, base.lat_name, base.time_name}
        self._check_variable_consistency(
            base.lon_name, base.lat_name, base.time_name
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

        meta = TopoInterrogator(path, var_name='z').interrogate_topo()
        with open('topo.data', 'a') as f:
            DescriptorWriter.write_topo_descriptor(f, meta)

    Usage — met (storm file)::

        meta = MetInterrogator(path, var_map).interrogate_met()
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
            Output of ``TopoInterrogator.interrogate_topo()``.
        """
        f.write(f"var_name       = {meta.var_name}\n")
        f.write(f"lon_name       = {meta.lon_name}\n")
        f.write(f"lat_name       = {meta.lat_name}\n")
        f.write(f"lon_convention = {meta.lon_convention}\n")
        f.write(f"lat_order      = {meta.lat_order}\n")
        f.write(f"dim_order      = {','.join(meta.dim_order)}\n")
        if meta.fill_value is not None:
            f.write(f"fill_value     = {meta.fill_value!r}\n")
        f.write(f"fill_action    = {meta.fill_action}\n")
        if meta.crop_bounds is not None:
            lon0, lon1, lat0, lat1 = meta.crop_bounds
            f.write(f"crop_bounds    = {lon0} {lon1} {lat0} {lat1}\n")
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
            Output of ``MetInterrogator.interrogate_met()``.
        """
        f.write("&file_info\n")
        f.write(f"  lon_name       = {meta.lon_name}\n")
        f.write(f"  lat_name       = {meta.lat_name}\n")
        f.write(f"  time_name      = {meta.time_name}\n")
        f.write(f"  dim_order      = {','.join(meta.dim_order)}\n")
        f.write(f"  lon_convention = {meta.lon_convention}\n")
        f.write(f"  lat_order      = {meta.lat_order}\n")
        if meta.fill_value is not None:
            f.write(f"  fill_value     = {meta.fill_value!r}\n")
        f.write(f"  fill_action    = {meta.fill_action}\n")
        f.write(f"  time_offset    = {meta.time_offset!r}\n")
        if meta.crop_bounds is not None:
            lon0, lon1, lat0, lat1 = meta.crop_bounds
            f.write(f"  crop_bounds    = {lon0} {lon1} {lat0} {lat1}\n")
        f.write("/\n")
        for var in meta.variables:
            f.write(
                f"&variable_info"
                f"  var_name={var.var_name}"
                f"  geoclaw_role={var.geoclaw_role}"
                f"  /\n"
            )
