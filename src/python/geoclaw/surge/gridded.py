#!/usr/bin/env python

r"""Gridded meteorological forcing.

``GriddedMetForcing`` describes file-backed wind/pressure fields (OWI/ASCII and
NetCDF).  It owns the file paths, the numeric ``file_format`` tag, the NetCDF
variable/dimension mapping and interpolation metadata, and the forcing controls
(scaling, storm time scaling, spatial crop, edge ramping, registration shift).

The descriptor reader (``read_data``) and writer (``write_data``) are moved
verbatim from the original ``Storm`` implementation so the emitted descriptor
bytes are byte-for-byte identical.  ``write_data`` reuses
``clawpack.geoclaw.netcdf_utils.MetInspector`` / ``DescriptorWriter``.
"""

from pathlib import Path

import numpy as np

from clawpack.geoclaw.surge.track import _Meta


class GriddedMetForcing(object):
    r"""Forcing read from external gridded wind/pressure field datasets."""

    # Units documented by each gridded storm format, used as a fallback for a
    # NetCDF met variable that carries no CF 'units' attribute.  OWI/NWS13
    # store pressure in millibar and wind in m/s (see the OWI header labels
    # "... Pressure Output in mb" / "... Wind Output ... in m/s").  Formats not
    # listed (netcdf/ERA5, which carry their own units) get no fallback.
    _met_format_units = {
        "nws13": {"wind_u": "m/s", "wind_v": "m/s", "pressure": "mbar"},
        "owi":   {"wind_u": "m/s", "wind_v": "m/s", "pressure": "mbar"},
    }

    def __init__(self, meta=None):
        self.meta = meta if meta is not None else _Meta()

        # Run-time modifications for storm
        self.scaling = [1.0, 1.0]          # Scaling of wind and pressure
        self.storm_time_scale = 1.0        # >1 slower storm, <1 faster storm

        # Spatial cropping + edge ramping - only applies to data storms.
        # crop_extent = [lon0, lon1, lat0, lat1] restricts the forcing to
        # a sub-region read directly from the file (NetCDF: strided read;
        # OWI/ASCII: read full then subset).  ramp_width tapers the wind and
        # pressure forcing to ambient over this many degrees inside the crop
        # (equivalently the file) edges.  crop_extent=None reads the full
        # file extent.
        self.crop_extent = None
        self.ramp_width = 1e0               # 1 degree

        # Registration shift (domain = file + shift) of the forcing grid,
        # parallel to topo/dtopo x_shift, y_shift.  Lets a projected or
        # mis-registered met grid be placed onto the domain.  0 = no shift.
        self.x_shift = 0.0
        self.y_shift = 0.0

        # NetCDF met-forcing metadata (populated by read_data for file_format==2)
        self.met_x_name = None
        self.met_y_name = None
        self.met_time_name = None
        self.met_lon_wrap = None
        self.met_y_increasing = None
        self.met_fill_value = None
        self.met_fill_action = None
        self.met_time_offset = None
        self.met_variable_map = {}

    # Shared bookkeeping (proxied through ``_Meta``).
    @property
    def time_offset(self):
        return self.meta.time_offset

    @time_offset.setter
    def time_offset(self, value):
        self.meta.time_offset = value

    @property
    def file_paths(self):
        return self.meta.file_paths

    @file_paths.setter
    def file_paths(self, value):
        self.meta.file_paths = value

    @property
    def file_format(self):
        return self.meta.file_format

    @file_format.setter
    def file_format(self, value):
        self.meta.file_format = value

    # =========================================================================
    # Read Routines
    @classmethod
    def read_data(cls, path: Path, verbose: bool=False):
        r"""Read in a data derived storm data information file

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
        """
        self = cls()

        with path.open() as data_file:
            data_file.readline()
            self.time_offset = np.datetime64(data_file.readline()[:19])
            self.file_format = int(
                                data_file.readline().partition("#")[0].rstrip())
            num_files = int(data_file.readline().partition("#")[0].rstrip())
            self.scaling = np.array([float(value)
                for value in data_file.readline().partition("#")[0].rstrip().split()],
                dtype=float)
            self.ramp_width = float(
                                data_file.readline().partition("#")[0].rstrip())
            _crop_vals = [float(v) for v in
                          data_file.readline().partition("#")[0].split()]
            # All-zero is the "no crop" sentinel (a valid extent requires
            # lon0<lon1 and lat0<lat1).
            if len(_crop_vals) == 4 and any(v != 0.0 for v in _crop_vals):
                self.crop_extent = _crop_vals
            else:
                self.crop_extent = None
            self.storm_time_scale = float(
                                data_file.readline().partition("#")[0].rstrip())
            self.x_shift = float(
                                data_file.readline().partition("#")[0].rstrip())
            self.y_shift = float(
                                data_file.readline().partition("#")[0].rstrip())
            data_file.readline()  # "# Format Data Information"

            if self.file_format == 1:
                # Skip blank line then "# File paths" comment
                data_file.readline()
                data_file.readline()
            elif self.file_format == 2:
                # Parse the &file_info / &variable_info descriptor produced
                # by DescriptorWriter.write_met_descriptor.
                file_info = {}
                variable_infos = []

                # Read &file_info block (lines until standalone "/")
                line = data_file.readline()
                if line.strip().startswith('&file_info'):
                    while True:
                        inner = data_file.readline()
                        if not inner:
                            break
                        s = inner.strip()
                        if s == '/':
                            break
                        if '=' in s and not s.startswith('&'):
                            key, _, val = s.partition('=')
                            file_info[key.strip()] = val.strip()

                # Read &variable_info lines until "# File paths" terminator
                while True:
                    line = data_file.readline()
                    if not line:
                        break
                    s = line.strip()
                    if s.startswith('# File paths'):
                        break
                    if s.startswith('&variable_info'):
                        var_info = {}
                        content = s[len('&variable_info'):].rstrip('/').strip()
                        for pair in content.split():
                            k, _, v = pair.partition('=')
                            if k and v:
                                var_info[k] = v
                        variable_infos.append(var_info)

                # Store fields on Storm object for Fortran to consume
                self.met_x_name = file_info.get('x_name')
                self.met_y_name = file_info.get('y_name')
                self.met_time_name = file_info.get('time_name')
                _conv = file_info.get('lon_wrap')
                self.met_lon_wrap = int(_conv) if _conv is not None else None
                _yi = file_info.get('y_increasing')
                self.met_y_increasing = (_yi == 'True') if _yi is not None else None
                _fv = file_info.get('fill_value')
                self.met_fill_value = float(_fv) if _fv is not None else None
                self.met_fill_action = file_info.get('fill_action', 'warn')
                _to = file_info.get('time_offset')
                self.met_time_offset = float(_to) if _to is not None else None
                self.met_variable_map = {
                    vi['geoclaw_role']: vi['var_name']
                    for vi in variable_infos
                    if 'geoclaw_role' in vi and 'var_name' in vi
                }
            else:
                raise TypeError(f"Unknown storm data file format type" +
                                f" '{self.file_format}' provided.")

            for i in range(num_files):
                self.file_paths.append(Path(data_file.readline().rstrip()))

        return self

    # =========================================================================
    # Write Routines
    def write_data(self, path, dim_mapping=None, var_mapping=None,
                   met_inspector=None, verbose=False):
        r"""
         """
        # Only one format right now
        _data_file_format_mapping = {'ascii': 1, 'nws12': 1, "owi": 1,
                                     'netcdf': 2, 'nws13': 2}

        if isinstance(self.file_format, int):
            file_format = self.file_format
        elif isinstance(self.file_format, str):
            if (self.file_format.lower() in
                                    _data_file_format_mapping.keys()):
                file_format = _data_file_format_mapping[
                                          self.file_format.lower()]
            else:
                raise TypeError(f"Unknown storm data file format type" +
                                f" '{self.file_format}' provided.")
        else:
            raise TypeError(f"Unknown storm data file format type" +
                            f" '{self.file_format}' provided.")

        # crop_extent: write a 4-value crop region or the all-zero "no
        # crop" sentinel (a valid extent requires lon0<lon1 and lat0<lat1).
        if self.crop_extent is None:
            crop_str = "0. 0. 0. 0."
        else:
            if len(self.crop_extent) != 4:
                raise ValueError(
                    "crop_extent must be [lon0, lon1, lat0, lat1], got "
                    f"{self.crop_extent!r}")
            crop_str = ' '.join(repr(float(v)) for v in self.crop_extent)

        with path.open("w") as data_file:
            # Write header
            data_file.write("# Data Derived Storm\n")

            # Time offset
            self.time_offset = np.datetime64(self.time_offset)
            if isinstance(self.time_offset, np.datetime64):
                t = np.datetime_as_string(self.time_offset, unit="s")
                data_file.write(f"{t.ljust(20)} # Time Offset\n")
            else:
                raise ValueError("Time offset must be a datetime64 object.")
            data_file.write(f"{str(file_format).ljust(20)} # File format\n")
            data_file.write(f"{str(len(self.file_paths)).ljust(20)} # Number of files\n")
            data_file.write(f"{str(self.scaling)[1:-1].replace(',', '').ljust(20)} # Scaling\n")
            data_file.write(f"{str(self.ramp_width).ljust(20)} # Ramp width\n")
            data_file.write(f"{crop_str.ljust(20)} # Met crop extent [lon0 lon1 lat0 lat1]; all-zero = full file\n")
            data_file.write(f"{str(self.storm_time_scale).ljust(20)} # Storm time scale (>1 slower, <1 faster)\n")
            data_file.write(f"{str(self.x_shift).ljust(20)} # Met x_shift (domain = file + shift)\n")
            data_file.write(f"{str(self.y_shift).ljust(20)} # Met y_shift (domain = file + shift)\n")
            data_file.write("# Format Data Information\n")
            if file_format == 1:
                # Check number of file paths
                if len(self.file_paths)%2 != 0:
                        raise ValueError("The number of files should be even, " +
                                         "one for pressure and wind, for each " +
                                         "resolution provided.")
                data_file.write("\n")
            elif file_format == 2:
                if len(self.file_paths) != 1:
                    raise ValueError(f"Expected 1 path for NetCDF format, " +
                                     f"got {len(self.file_paths)}")

                from clawpack.geoclaw.netcdf_utils import (
                    MetInspector, DescriptorWriter)

                if met_inspector is None:
                    # MetInspector auto-discovers met roles (CF standard_name
                    # first, then common variable names); explicit var_mapping
                    # entries override individual roles, so callers only need
                    # to supply the non-standard names.  Coordinates are
                    # discovered via CF axis/standard_name conventions;
                    # dim_mapping is accepted for backwards compatibility but
                    # not needed.  format_units supplies the storm format's
                    # documented units (e.g. NWS13/OWI pressure = mbar) for any
                    # variable that carries no CF 'units' attribute.
                    _fmt_units = self._met_format_units.get(
                        str(self.file_format).lower())
                    with MetInspector(self.file_paths[0],
                                         variable_map=var_mapping,
                                         time_reference=self.time_offset,
                                         format_units=_fmt_units) as mi:
                        meta = mi.inspect_met()
                else:
                    meta = met_inspector.inspect_met()

                # The spatial crop is carried by the top-level
                # crop_extent line (read by set_storm and applied at file
                # read for both formats), so it is intentionally not also
                # emitted in the NetCDF descriptor's crop_bounds.
                DescriptorWriter.write_met_descriptor(data_file, meta)

            # Write paths
            data_file.write("# File paths\n")
            for path in self.file_paths:
                data_file.write(f"{path}\n")
