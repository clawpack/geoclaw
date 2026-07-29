#!/usr/bin/env python

r"""Oceanweather (OWI) WIN/PRE meteorological-forcing field I/O.

This module reads and writes the fixed-width Oceanweather ASCII wind/pressure
format (a.k.a. NWS12/OWI) used by ADCIRC and USACE.  It is the Python-side
counterpart to the Fortran reader ``read_OWI_ASCII`` /
``read_OWI_ASCII_header`` in ``gridded_met_forcing_module.f90`` and produces
files that reader consumes.

Format contract (80-column fixed width)
---------------------------------------
Each of the pressure (``.PRE``) and wind (``.WIN``) files begins with a single
file header line::

    Oceanweather WIN/PRE Format<pad>YYYYMMDDHH     YYYYMMDDHH

where the start date begins in column 56 and the end date in column 71 (the
Fortran header read is ``(t56,i4,i2,i2,i2,t71,i4,i2,i2,i2)``).

Each time step is then a grid-spec header line followed by the field data::

    iLat=<i4>iLong=<i4>DX=<f6.4>DY=<f6.4>SWLat=<f8.4>SWLon=<f8.4>DT=YYYYMMDDHHMM

with ``iLat`` = number of latitude points, ``iLong`` = number of longitude
points, and ``DT`` at column 69.  Field values follow, eight per line in
``f10.4``, ordered longitude-fastest (``((f(i,j), i=1..iLong), j=1..iLat)``).
The ``.PRE`` file stores one pressure block per time; the ``.WIN`` file stores
the ``u`` block then the ``v`` block per time.

Units
-----
Pressure is stored in the file in millibar (the Fortran reader multiplies by
100 to obtain Pa).  This module's in-memory representation is SI: ``read_owi``
returns pressure in **Pa** and ``write_owi`` expects Pa, handling the mb<->Pa
conversion at the file boundary.  Wind is m/s in both the file and memory.

Precision: field values are written as fixed-width ``f10.4``, so a write/read
round-trip preserves pressure only to ~1e-2 Pa (1e-4 mbar) and wind to
~1e-4 m/s -- ample for forcing, but not bit-exact.  (Tests assert the pressure
round-trip to ``atol=1e-2`` Pa for this reason.)

Coordinates
-----------
The SW corner (``SWLon``, ``SWLat``) is the first grid point, so
``longitude = SWLon + arange(iLong) * DX`` and likewise for latitude.  (Note:
the Fortran reader uses a one-based ``sw + i*dx`` that omits the SW corner --
a pre-existing off-by-one there; this reader follows the OWI convention.)
"""

from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np

# File header layout: label left-justified into this width, then the start and
# end dates separated by five spaces, so the start date lands in column 56.
_HEADER_LABEL = "Oceanweather WIN/PRE Format"
_HEADER_LABEL_WIDTH = 55
_FILE_DATE_FMT = "%Y%m%d%H"        # file header start/end (hour precision)
_RECORD_DATE_FMT = "%Y%m%d%H%M"    # per-record DT (minute precision)
_VALUES_PER_LINE = 8
_PA_PER_MBAR = 100.0


@dataclass
class OWIData:
    r"""In-memory representation of an OWI WIN/PRE dataset.

    :Attributes:
     - *time* (numpy.ndarray) ``datetime64[s]`` array of length ``nt``.
     - *longitude* (numpy.ndarray) 1D longitudes (length ``nx``), ascending.
     - *latitude* (numpy.ndarray) 1D latitudes (length ``ny``), ascending.
     - *wind_u*, *wind_v* (numpy.ndarray) ``(nt, ny, nx)`` wind components (m/s).
     - *pressure* (numpy.ndarray) ``(nt, ny, nx)`` pressure (Pa).
    """
    time: np.ndarray
    longitude: np.ndarray
    latitude: np.ndarray
    wind_u: np.ndarray
    wind_v: np.ndarray
    pressure: np.ndarray

    @property
    def num_times(self):
        return self.time.shape[0]

    @property
    def shape(self):
        r"""Spatial grid shape ``(ny, nx)``."""
        return (self.latitude.shape[0], self.longitude.shape[0])


def _as_datetime(value):
    r"""Coerce a datetime64 / datetime / string to a naive ``datetime``."""
    dt64 = np.datetime64(value, "s")
    # datetime64[s] epoch is 1970-01-01T00:00:00; build an aware UTC datetime
    # then drop tzinfo for strftime (the format carries no zone).
    epoch = np.datetime64("1970-01-01T00:00:00", "s")
    seconds = int((dt64 - epoch) / np.timedelta64(1, "s"))
    return datetime.fromtimestamp(seconds, tz=timezone.utc).replace(tzinfo=None)


def _parse_record_header(line):
    r"""Parse a grid-spec header line into (ny, nx, dx, dy, swlat, swlon, dt).

    Robust to the small formatting variations seen in real OWI files (e.g.
    ``SWLat`` written with 4 or 5 decimals) by splitting on the ``key=value``
    tokens rather than fixed columns.
    """
    def _field(key, width):
        start = line.index(key) + len(key)
        return line[start:start + width]

    ny = int(_field("iLat=", 4))
    nx = int(_field("iLong=", 4))
    dx = float(_field("DX=", 6))
    dy = float(_field("DY=", 6))
    swlat = float(_field("SWLat=", 8))
    swlon = float(_field("SWLon=", 8))
    dt = datetime.strptime(_field("DT=", 12), _RECORD_DATE_FMT)
    return ny, nx, dx, dy, swlat, swlon, dt


def _read_blocks(fh, num_blocks, count):
    r"""Read ``num_blocks`` flat field blocks of ``count`` values each.

    Each block is ``ceil(count / 8)`` lines of up to eight ``f10.4`` values.
    Returns a list of ``num_blocks`` 1D float arrays.
    """
    blocks = []
    for _ in range(num_blocks):
        values = []
        while len(values) < count:
            line = fh.readline()
            if not line:
                raise EOFError("Unexpected end of OWI file while reading data.")
            values.extend(float(v) for v in line.split())
        blocks.append(np.array(values[:count], dtype=float))
    return blocks


def read_owi(pressure_path, wind_path):
    r"""Read an OWI WIN/PRE pair into an :class:`OWIData` object.

    :Input:
     - *pressure_path* (path-like) the ``.PRE`` file.
     - *wind_path* (path-like) the ``.WIN`` file.

    :Output:
     - (:class:`OWIData`) fields with pressure in Pa and wind in m/s, arrays
       shaped ``(nt, ny, nx)``.
    """
    times = []
    u_fields, v_fields, p_fields = [], [], []
    lon = lat = None

    with open(pressure_path, "r") as pre, open(wind_path, "r") as win:
        pre.readline()   # file header (start/end dates -- not needed here)
        win.readline()

        while True:
            win_header = win.readline()
            if not win_header or not win_header.strip():
                break
            pre.readline()   # pressure record header (dims taken from wind)

            ny, nx, dx, dy, swlat, swlon, dt = _parse_record_header(win_header)
            count = nx * ny

            if lon is None:
                lon = swlon + np.arange(nx) * dx
                lat = swlat + np.arange(ny) * dy

            u_flat, v_flat = _read_blocks(win, 2, count)
            (p_flat,) = _read_blocks(pre, 1, count)

            times.append(np.datetime64(dt, "s"))
            u_fields.append(u_flat.reshape(ny, nx))
            v_fields.append(v_flat.reshape(ny, nx))
            # File stores millibar; convert to Pa for the SI in-memory form.
            p_fields.append(p_flat.reshape(ny, nx) * _PA_PER_MBAR)

    return OWIData(time=np.array(times, dtype="datetime64[s]"),
                   longitude=lon, latitude=lat,
                   wind_u=np.array(u_fields),
                   wind_v=np.array(v_fields),
                   pressure=np.array(p_fields))


def read_owi_start_time(path):
    r"""Return the first record's timestamp as ``datetime64[s]``.

    Cheap: reads only the file header and first grid-spec header, so it is
    suitable for inferring a descriptor ``time_offset`` without loading the
    (potentially large) field data.
    """
    with open(path, "r") as fh:
        fh.readline()                       # file header
        _, _, _, _, _, _, dt = _parse_record_header(fh.readline())
    return np.datetime64(dt, "s")


def _format_record_header(ny, nx, dx, dy, swlat, swlon, dt):
    r"""Build one 80-column grid-spec header line (see module docstring)."""
    return (f"iLat={ny:4d}iLong={nx:4d}DX={dx:6.4f}DY={dy:6.4f}"
            f"SWLat={swlat:8.4f}SWLon={swlon:8.4f}"
            f"DT={dt.strftime(_RECORD_DATE_FMT)}")


def _write_blocks(fh, *flats):
    r"""Write flat field arrays as eight ``f10.4`` values per line."""
    for flat in flats:
        for i, value in enumerate(flat):
            fh.write(f"{value:10.4f}")
            if (i + 1) % _VALUES_PER_LINE == 0:
                fh.write("\n")
        if len(flat) % _VALUES_PER_LINE != 0:
            fh.write("\n")


def write_owi(data, pressure_path, wind_path):
    r"""Write an :class:`OWIData` object to an OWI WIN/PRE pair.

    Produces files readable by the Fortran ``read_OWI_ASCII`` path: pressure is
    converted from Pa back to millibar, wind is written in m/s, and values are
    ordered longitude-fastest.

    :Input:
     - *data* (:class:`OWIData`) fields to write (Pa / m/s, ``(nt, ny, nx)``).
     - *pressure_path* (path-like) output ``.PRE`` file.
     - *wind_path* (path-like) output ``.WIN`` file.
    """
    ny, nx = data.shape
    dx = float(data.longitude[1] - data.longitude[0]) if nx > 1 else 0.0
    dy = float(data.latitude[1] - data.latitude[0]) if ny > 1 else 0.0
    swlon = float(data.longitude[0])
    swlat = float(data.latitude[0])

    start = _as_datetime(data.time[0])
    end = _as_datetime(data.time[-1])
    file_header = (f"{_HEADER_LABEL.ljust(_HEADER_LABEL_WIDTH)}"
                   f"{start.strftime(_FILE_DATE_FMT)}     "
                   f"{end.strftime(_FILE_DATE_FMT)}")

    with open(pressure_path, "w") as pre, open(wind_path, "w") as win:
        pre.write(file_header + "\n")
        win.write(file_header + "\n")

        for n in range(data.num_times):
            dt = _as_datetime(data.time[n])
            header = _format_record_header(ny, nx, dx, dy, swlat, swlon, dt)
            pre.write(header + "\n")
            win.write(header + "\n")

            _write_blocks(win, data.wind_u[n].ravel(), data.wind_v[n].ravel())
            # Convert Pa back to the millibar the file format stores.
            _write_blocks(pre, (data.pressure[n] / _PA_PER_MBAR).ravel())
