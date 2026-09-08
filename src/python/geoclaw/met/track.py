#!/usr/bin/env python

r"""Track and StormTrack objects for meteorological forcing.

This module defines the ``Track`` and ``StormTrack`` objects used by the
meteorological-forcing refactor along with the track-format readers that
produce a ``StormTrack`` (ATCF, HURDAT, IBTrACS, JMA, TCVITALS).  The parsing
logic here is moved verbatim from the original ``Storm`` readers so the parsed
values (and any downstream emitted bytes) are byte-for-byte identical.

A ``Track`` is a generic evolving spatial feature (a center location over
time).  A ``StormTrack`` adds storm-specific metadata (intensity, radii,
classification, basin).  The time axis is a NumPy ``datetime64`` array,
uniform with the gridded field axis (see ``met_forcing_refactor.md`` Section
6).
"""

import datetime
import enum
import re
import warnings

import numpy as np
import pandas as pd

import clawpack.geoclaw.units as units


# =============================================================================
#  Common acronyms across formats

# ATCF basins with their expanded names
# see https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
ATCF_basins = {"AL": "Atlantic",
               "CP": "Central Pacific",
               "EP": "East Pacific",
               "IO": "North Indian Ocean",
               "SH": "Southern Hemisphere",
               "SL": "Southern Atlantic",
               "LS": "Southern Atlantic",
               "WP": "North West Pacific"}

# TCVitals basins with their expanded names
# see http://www.emc.ncep.noaa.gov/HWRF/tcvitals-draft.html
TCVitals_Basins = {"L": "North Atlantic",
                   "E": "North East Pacific",
                   "C": "North Central Pacific",
                   "W": "North West Pacific",
                   "B": "Bay of Bengal (North Indian Ocean)",
                   "A": "Arabian Sea (North Indian Ocean)",
                   "Q": "South Atlantic",
                   "P": "South Pacific",
                   "S": "South Indian Ocean"}

# Tropical Cyclone Designations
# see https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
TC_designations = {"DB": "disturbance",
                   "TD": "tropical depression",
                   "TS": "tropical storm",
                   "TY": "typhoon",
                   "ST": "super typhoon",
                   "TC": "tropical cyclone",
                   "HU": "hurricane",
                   "SD": "subtropical depression",
                   "SS": "subtropical storm",
                   "EX": "extratropical systems",
                   "IN": "inland",
                   "DS": "dissipating",
                   "LO": "low",
                   "WV": "tropical wave",
                   "ET": "extrapolated",
                   "XX": "unknown"}

# HURDAT special designations
# see http://www.aoml.noaa.gov/hrd/data_sub/newHURDAT.html
hurdat_special_entries = {"L": "landfall",
                          "W": "max wind",
                          "P": "min pressure",
                          "I": "max intensity",
                          "C": "closest approach",
                          "S": "status change",
                          "G": "genesis",
                          "T": "additional track point"}

# ---------------------------------------------------------------------------
# Standardized vocabularies
#
# Each archive spells basin and storm status its own way, and the maps above
# expand them into prose that differs between formats ("Atlantic" from ATCF vs
# "North Atlantic" from tcvitals).  These enums give one canonical value per
# concept so code can compare across formats.
#
# The string attributes ``StormTrack.basin`` and ``.classification`` are
# unchanged -- these are additive, reached through ``basin_code`` and
# ``status()``.  Turning the existing attributes into enums would change what
# every reader reports and how it serializes, for no gain to code that already
# works.
# ---------------------------------------------------------------------------

class Basin(enum.Enum):
    r"""Canonical tropical-cyclone basin, independent of source vocabulary."""
    NORTH_ATLANTIC = "north_atlantic"
    SOUTH_ATLANTIC = "south_atlantic"
    EAST_PACIFIC = "east_pacific"
    CENTRAL_PACIFIC = "central_pacific"
    WEST_PACIFIC = "west_pacific"
    NORTH_INDIAN = "north_indian"
    SOUTH_INDIAN = "south_indian"
    SOUTH_PACIFIC = "south_pacific"
    SOUTHERN_HEMISPHERE = "southern_hemisphere"
    UNKNOWN = "unknown"


class StormStatus(enum.Enum):
    r"""Canonical storm status (development stage), not Saffir-Simpson category."""
    DISTURBANCE = "disturbance"
    TROPICAL_WAVE = "tropical_wave"
    LOW = "low"
    TROPICAL_DEPRESSION = "tropical_depression"
    TROPICAL_STORM = "tropical_storm"
    HURRICANE = "hurricane"
    TYPHOON = "typhoon"
    SUPER_TYPHOON = "super_typhoon"
    TROPICAL_CYCLONE = "tropical_cyclone"
    SUBTROPICAL_DEPRESSION = "subtropical_depression"
    SUBTROPICAL_STORM = "subtropical_storm"
    EXTRATROPICAL = "extratropical"
    DISSIPATING = "dissipating"
    INLAND = "inland"
    UNKNOWN = "unknown"


# Source code -> Basin, per format.  ATCF/HURDAT2 share the two-letter codes;
# IBTrACS uses its own; tcvitals uses single letters.
_BASIN_CODES = {
    "atcf": {"AL": Basin.NORTH_ATLANTIC, "CP": Basin.CENTRAL_PACIFIC,
             "EP": Basin.EAST_PACIFIC, "IO": Basin.NORTH_INDIAN,
             "SH": Basin.SOUTHERN_HEMISPHERE, "SL": Basin.SOUTH_ATLANTIC,
             "LS": Basin.SOUTH_ATLANTIC, "WP": Basin.WEST_PACIFIC},
    "ibtracs": {"NA": Basin.NORTH_ATLANTIC, "SA": Basin.SOUTH_ATLANTIC,
                "EP": Basin.EAST_PACIFIC, "WP": Basin.WEST_PACIFIC,
                "NI": Basin.NORTH_INDIAN, "SI": Basin.SOUTH_INDIAN,
                "SP": Basin.SOUTH_PACIFIC},
    "tcvitals": {"L": Basin.NORTH_ATLANTIC, "E": Basin.EAST_PACIFIC,
                 "C": Basin.CENTRAL_PACIFIC, "W": Basin.WEST_PACIFIC,
                 "B": Basin.NORTH_INDIAN, "A": Basin.NORTH_INDIAN,
                 "Q": Basin.SOUTH_ATLANTIC, "P": Basin.SOUTH_PACIFIC,
                 "S": Basin.SOUTH_INDIAN},
}
_BASIN_CODES["hurdat"] = _BASIN_CODES["atcf"]

# Source code -> StormStatus.  ATCF/HURDAT2/IBTrACS all use the TC designation
# codes, with IBTrACS adding a few of its own.
_STATUS_CODES = {
    "DB": StormStatus.DISTURBANCE, "WV": StormStatus.TROPICAL_WAVE,
    "LO": StormStatus.LOW, "TD": StormStatus.TROPICAL_DEPRESSION,
    "TS": StormStatus.TROPICAL_STORM, "HU": StormStatus.HURRICANE,
    "HR": StormStatus.HURRICANE, "TY": StormStatus.TYPHOON,
    "ST": StormStatus.SUPER_TYPHOON, "TC": StormStatus.TROPICAL_CYCLONE,
    "SD": StormStatus.SUBTROPICAL_DEPRESSION,
    "SS": StormStatus.SUBTROPICAL_STORM, "SU": StormStatus.SUBTROPICAL_STORM,
    "EX": StormStatus.EXTRATROPICAL, "ET": StormStatus.EXTRATROPICAL,
    "DS": StormStatus.DISSIPATING, "IN": StormStatus.INLAND,
    "MD": StormStatus.LOW, "PT": StormStatus.LOW,
    "XX": StormStatus.UNKNOWN, "NR": StormStatus.UNKNOWN,
}


def standardize_basin(code, source="atcf"):
    r"""Map a format's basin *code* onto a :class:`Basin`.

    Unrecognized codes give ``Basin.UNKNOWN`` rather than raising: a new or
    misspelled code in an archive should not stop a track from being read.
    """
    if code is None:
        return Basin.UNKNOWN
    table = _BASIN_CODES.get(str(source).lower(), {})
    return table.get(str(code).strip().upper(), Basin.UNKNOWN)


def standardize_status(code, source="atcf"):
    r"""Map a format's status *code* onto a :class:`StormStatus`.

    Unrecognized codes give ``StormStatus.UNKNOWN`` rather than raising.
    """
    if code is None:
        return StormStatus.UNKNOWN
    return _STATUS_CODES.get(str(code).strip().upper(), StormStatus.UNKNOWN)


# Quadrant order used by ATCF, HURDAT2 and IBTrACS alike.
WIND_RADII_QUADRANTS = ("NE", "SE", "SW", "NW")

# Wind-speed thresholds (knots) at which quadrant radii are reported.
WIND_RADII_THRESHOLDS_KT = (34.0, 50.0, 64.0)


# Warning for formats that have yet to have a default way to determine crticial
# radii from the input data
missing_data_warning_str = """*** Cannot yet automatically determine the
    maximum wind radius.  Will write out GeoClaw
    formats but note that these will not work
    when running GeoClaw currently without a custom
    `max_wind_radius_fill` function passed as argument
    to the `write` function."""

# Warning for not having any time points with both a max wind speed and central
# pressure observation
missing_necessary_data_warning_str = """No storm points in the input file
    had both a max wind speed and a central pressure observation."""


class NoDataError(ValueError):
    """Exception to raise when no valid data in input file"""
    pass


# ---------------------------------------------------------------------------
# Missing-data contract
#
# In memory, ``np.nan`` is the *only* marker for a missing value, in every reader
# and every field.  Archive-specific sentinels are normalized here, at the reader
# boundary, so nothing downstream has to know which format a track came from.
# ``read_atcf`` has always worked this way; hurdat/jma/ibtracs used ``-1``, which
# silently defeated ``write_geoclaw``'s NaN-based fill and skip logic.
#
# Note that zero is *not* missing.  HURDAT2 reports genuinely-zero quadrant wind
# radii for weak systems (no 34-kt winds exist to have a radius), which is real
# data and distinct from its ``-999`` "unknown".
# ---------------------------------------------------------------------------

# AOML HURDAT2 uses -99 for unknown wind and -999 for unknown pressure/radii.
HURDAT_SENTINELS = (-99.0, -999.0, -9999.0)

# NCEP tcvitals uses -99/-999 in the same role.
TCVITALS_SENTINELS = (-99.0, -999.0)


# Default order in which agencies are consulted for IBTrACS intensity fields
# when the `wmo_*` variable is missing and neither `wmo_agency` nor `usa_agency`
# names a provider.  Shared by StormTrack.read_ibtracs and iter_ibtracs so the
# two cannot drift apart.
_IBTRACS_AGENCY_PREF = ('wmo', 'usa', 'tokyo', 'newdelhi', 'reunion', 'bom',
                        'nadi', 'wellington', 'cma', 'hko', 'ds824', 'td9636',
                        'td9635', 'neumann', 'mlc')


def _atcf_quadrant_radii(df):
    r"""Quadrant wind radii from a raw ATCF frame, indexed by ``DATE``.

    ATCF records each wind threshold separately at the same ``DATE``, so the
    thresholds must be pivoted out *before* the reader collapses each date to a
    single row.  Returns a DataFrame indexed by ``DATE`` with 12 columns ordered
    ``(threshold, quadrant)`` -- 34/50/64 kt x NE/SE/SW/NW -- in nautical miles.

    A zero is treated as missing here, matching how this reader already handles
    ATCF's other radius fields: the format leaves an unavailable quadrant as 0,
    so a 0 cannot be distinguished from a genuine zero radius.  (HURDAT2 *can*
    distinguish them -- it has an explicit -999 -- so there a zero is kept.)
    """
    columns = ["RAD1", "RAD2", "RAD3", "RAD4"]
    frame = df[["DATE", "RAD"] + columns].copy()
    frame[columns] = frame[columns].where(frame[columns] != 0, np.nan)

    pivot = frame.pivot_table(index="DATE", columns="RAD", values=columns,
                              aggfunc="first")
    # Reorder to (threshold, quadrant) and fill in any threshold the storm never
    # reached, so the shape is the same for every track.
    wanted = [(column, threshold)
              for threshold in WIND_RADII_THRESHOLDS_KT
              for column in columns]
    return pivot.reindex(columns=pd.MultiIndex.from_tuples(wanted))


def _sentinel_to_nan(value, sentinels):
    r"""Replace archive missing-value *sentinels* with ``np.nan``.

    Works on a scalar or an array.  Applied *before* any unit conversion, so a
    sentinel is never scaled into a plausible-looking physical value -- e.g.
    ``units.convert(-999, 'mbar', 'Pa')`` would otherwise yield -99900.0 Pa.

    An empty or whitespace-only field is also missing: fixed-width archives pad
    absent values with blanks, and a HURDAT2 revision that predates a column
    simply ends the line early (the bundled fixture stops at 20 fields, where
    current files carry 21).
    """
    values = value if isinstance(value, (list, tuple)) else [value]
    cleaned = []
    for item in values:
        if isinstance(item, str) and not item.strip():
            cleaned.append(np.nan)
        else:
            cleaned.append(item)

    array = np.asarray(cleaned, dtype=float)
    missing = np.zeros(array.shape, dtype=bool)
    for sentinel in sentinels:
        missing |= (array == sentinel)
    result = np.where(missing, np.nan, array)
    if not isinstance(value, (list, tuple)):
        return float(result[0])
    return result


class _Meta(object):
    r"""Shared meteorological-forcing bookkeeping.

    Holds the fields that are common to the track/parameterized path and the
    gridded path and that the readers populate alongside the track/dataset
    data: the Python-to-Fortran ``time_offset``, the list of ``file_paths``,
    and the ``file_format`` tag.  A single instance is shared between a
    ``StormTrack`` and the forcing objects that reference it so the values
    stay consistent regardless of which path is active.
    """

    def __init__(self):
        # Time offset (usually landfall) - see StormTrack/Storm docs.
        self.time_offset = None
        # Original file(s) read in (track) or pointed to (gridded).
        self.file_paths = []
        # Format read in or the type of data-driven storm.
        self.file_format = None


class Track(object):
    r"""Generic evolving feature with a spatial center over time.

    :Attributes:
     - *t* (ndarray) datetime64 array of times for each track point.
     - *center* (ndarray(:, :)) ``(n, 2)`` array of ``(lon, lat)``.  For a
       tropical cyclone this is the eye; for other systems it is the central
       point (e.g. a pressure minimum).
     - *name*, *id*, *event* optional metadata.
    """

    def __init__(self, t=None, center=None, name=None, id=None, event=None,
                 meta=None):
        self.t = t
        self.center = center
        self.name = name
        self.id = id
        self.event = event
        self.meta = meta if meta is not None else _Meta()

    # ``eye_location`` is the historical name for the (lon, lat) center of a
    # storm; keep it as an alias of ``center`` so verbatim reader/writer code
    # (and the legacy ``Storm`` surface) keep working.
    @property
    def eye_location(self):
        return self.center

    @eye_location.setter
    def eye_location(self, value):
        self.center = value

    # ``ID`` is the historical (upper-case) name of the track identifier.
    @property
    def ID(self):
        return self.id

    @ID.setter
    def ID(self, value):
        self.id = value

    # Shared bookkeeping proxied through ``meta`` so readers can set these on
    # ``self`` verbatim.
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


class StormTrack(Track):
    r"""A :class:`Track` carrying storm-specific metadata.

    Adds intensity-like parameters (``max_wind_speed``, ``max_wind_radius``,
    ``central_pressure``, ``storm_radius``) and optional descriptors
    (``classification``, ``basin``, ``wind_speeds``).  Supports both tropical
    and extratropical systems.  The track-format readers below are classmethods
    that produce a populated ``StormTrack``.

    Quadrant wind radii, where the source provides them, are in
    ``wind_radii`` with shape ``(n_times, n_thresholds, 4)`` in metres, the last
    axis ordered NE, SE, SW, NW (:data:`WIND_RADII_QUADRANTS`) and the middle
    axis matching ``wind_radii_thresholds`` (34/50/64 kt as m/s).  ``np.nan``
    means *not reported*.  Note the archives differ on whether a zero is data:
    HURDAT2 has an explicit ``-999`` for unknown, so a zero radius there is a
    real zero and is kept; ATCF leaves an unavailable quadrant as ``0``, which
    is indistinguishable from a genuine zero, so those become ``np.nan``.

    ``basin`` and ``classification`` remain the source's own strings.  For
    comparisons across formats use :meth:`basin_standard` and :meth:`status`,
    which map onto the :class:`Basin` and :class:`StormStatus` enums.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.max_wind_speed = None
        self.max_wind_radius = None
        self.central_pressure = None
        self.storm_radius = None
        self.classification = None
        self.basin = None
        self.basin_code = None
        self.wind_speeds = None
        self.wind_radii = None
        self.wind_radii_thresholds = None

    def basin_standard(self):
        r"""This track's basin as a :class:`Basin`, or ``Basin.UNKNOWN``."""
        return standardize_basin(self.basin_code, source=self.file_format
                                 or "atcf")

    def status(self):
        r"""Per-time storm status as :class:`StormStatus` values.

        ``None`` when the source carried no classification.
        """
        if self.classification is None:
            return None
        return np.array([standardize_status(code) for code
                         in np.atleast_1d(self.classification)], dtype=object)

    # =========================================================================
    # Read Routines (moved verbatim from Storm)
    @classmethod
    def read_atcf(cls, path, verbose=False):
        r"""Read in a ATCF formatted storm file

        ATCF format has storm stored individually so there is no support for
        multiple storms in a particular file.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
        """
        self = cls()

        # See here for the ATCF format documentation:
        #   https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abdeck.txt

        # Slightly more robust converter for ATCF data fields that can be
        # missing
        def num_converter(x):
            if isinstance(x, str):
                if len(x.strip()) == 0:
                    # Only whitespace
                    return np.nan
                else:
                    # Assume this is still a number
                    return float(x)
            elif x is None:
                return np.nan
            return float(x)

        df = pd.read_csv(path, engine="python", sep=",+", names=[
            "BASIN", "CY", "YYYYMMDDHH", "TECHNUM", "TECH", "TAU",
            "LAT", "LON", "VMAX", "MSLP", "TY",
            "RAD", "WINDCODE", "RAD1", "RAD2", "RAD3", "RAD4",
            "POUTER", "ROUTER", "RMW", "GUSTS", "EYE", "SUBREGION",
            "MAXSEAS", "INITIALS", "DIR", "SPEED", "STORMNAME", "DEPTH",
            "SEAS", "SEASCODE", "SEAS1", "SEAS2", "SEAS3", "SEAS4",
            "USERDEFINE1", "userdata1",
            "USERDEFINE2", "userdata2",
            "USERDEFINE3", "userdata3",
            "USERDEFINE4", "userdata4",
            "USERDEFINE5", "userdata5",
        ],
            converters={
                "YYYYMMDDHH": lambda d: np.datetime64(
                        f"{d[1:5]}-{d[5:7]}-{d[7:9]}T{d[9:11]}"),
                "TAU": lambda d: datetime.timedelta(hours=int(d)),
                "LAT": lambda d: (-.1 if d[-1] == "S" else .1) * int(d.strip("NS ")),
                "LON": lambda d: (-.1 if d[-1] == "W" else .1) * int(d.strip("WE ")),
                "RAD": num_converter,
                "RAD": num_converter,
                "RAD1": num_converter,
                "RAD2": num_converter,
                "RAD3": num_converter,
                "RAD4": num_converter,
                "ROUTER": num_converter,
                "RMW": num_converter,
                "STORMNAME": lambda d: (d.strip() if isinstance(d, str) else d)
        },
            dtype={
                "BASIN": str,
                "CY": int,
                "VMAX": float,
                "MSLP": float,
                "TY": str
        })

        # Grab data regarding basin and cyclone number from first row
        df_basin_code = df["BASIN"][0]
        self.basin = ATCF_basins[df_basin_code]
        self.ID = df["CY"][0]

        # Keep around the name as an array
        self.name = df["STORMNAME"].to_numpy()

        # Take forecast period TAU into consideration
        df['DATE'] = df["YYYYMMDDHH"] + df["TAU"]

        # Quadrant wind radii, built from the *pre-groupby* frame.  ATCF stores
        # each wind threshold as its own record at the same DATE, so the
        # groupby("DATE").first() below keeps only the 34-kt row and discards
        # the 50- and 64-kt ones.  This pivot is deliberately a separate pass so
        # the legacy pipeline below is untouched -- the committed
        # atcf_geoclaw.txt golden depends on it byte for byte.
        quadrant_radii = _atcf_quadrant_radii(df)

        df = df[["DATE", "TAU", "TY", "LAT", "LON", "VMAX", "MSLP",
                 "ROUTER", "RMW", "RAD", "RAD1", "RAD2", "RAD3", "RAD4", ]]
        df = df.sort_values(by=["DATE", "TAU"]).reset_index(drop=True)

        # For each DATE, choose best (smallest TAU) available data
        for c in ["LAT", "LON", "VMAX", "MSLP", "ROUTER", "RMW",
                  "RAD", "RAD1", "RAD2", "RAD3", "RAD4"]:
            df[c] = df[c].where(df[c] != 0, np.nan)  # value 0 means NaN
            df[c] = df.groupby("DATE")[c].bfill()
        df = df.groupby("DATE").first()

        # Wind profile (occasionally missing for older ATCF storms)
        # Wind speeds and their radii
        df["RAD_MEAN"] = df[["RAD1", "RAD2", "RAD3", "RAD4"]].mean(
            axis=1, skipna=True)
        df = df.drop(["TAU", "RAD1", "RAD2", "RAD3", "RAD4"], axis=1)
        df = df.dropna(how="any", subset=["LAT", "LON"])

        # Create time
        # self.t = list(df.index.to_pydatetime())
        self.t = df.index

        # Classification, note that this is not the category of the storm
        self.classification = df["TY"].to_numpy()

        # Eye location - longitude/latitude order
        self.eye_location = df[["LON", "LAT"]].to_numpy()

        # Convert to correct units:
        #  max_wind_speed - Convert knots to m/s - 0.51444444
        #  max_wind_radius  - convert from nm to m - 1.8520000031807990 * 1000.0
        #  central_pressure - convert from mbar to Pa - 100.0
        #  Radius of last isobar contour - convert from nm to m - 1.852000003180799d0 * 1000.0
        self.max_wind_speed = units.convert(
            df["VMAX"].to_numpy(), 'knots', 'm/s')
        self.central_pressure = units.convert(
            df["MSLP"].to_numpy(), 'mbar', 'Pa')
        self.max_wind_radius = units.convert(df["RMW"].to_numpy(), 'nmi', 'm')
        self.storm_radius = units.convert(df["ROUTER"].to_numpy(), 'nmi', 'm')
        self.wind_speeds = df[["RAD", "RAD_MEAN"]].to_numpy()
        self.wind_speeds[:, 0] = units.convert(
            self.wind_speeds[:, 0], 'knots', 'm/s')
        self.wind_speeds[:, 1] = units.convert(
            self.wind_speeds[:, 1], 'nmi', 'm')

        # Align the quadrant radii onto the surviving times.
        self.wind_radii = quadrant_radii.reindex(df.index).to_numpy().reshape(
            len(df.index), len(WIND_RADII_THRESHOLDS_KT), 4)
        self.wind_radii = units.convert(self.wind_radii, 'nmi', 'm')
        self.wind_radii_thresholds = units.convert(
            np.array(WIND_RADII_THRESHOLDS_KT), 'knots', 'm/s')

        self.basin_code = str(df_basin_code)

        self.file_paths.append(path)
        self.file_format = "atcf"

        return self

    @classmethod
    def read_hurdat(cls, path, verbose=False, storm_id=None, name=None,
                    year=None):
        r"""Read a HURDAT2 formatted storm file.

        This is the current version of HURDAT data available (HURDAT 2).  The
        archives NOAA/AOML distributes hold *every* storm in a basin in one file
        (~2000 storms for the Atlantic, 1851-present), so pass a selector to pick
        one; a file containing a single storm needs no selector.  Use
        :func:`iter_hurdat` to walk every storm, or :func:`catalog_hurdat` for an
        index of what a file contains.

        For more details on the HURDAT format and getting data see

        http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
         - *storm_id* (string) ATCF-style id to select, e.g. ``"AL092008"``.
         - *name* (string) Storm name to select, e.g. ``"IKE"``
           (case-insensitive).
         - *year* (int) Season to select; combine with *name* when a name has
           been reused across seasons.

        :Raises:
         - *ValueError* If no storm matches the selector, or if the file holds
           more than one storm and no selector narrows it to exactly one.
        """
        matches = [(header, lines) for (header, lines) in _hurdat_blocks(path)
                   if _hurdat_selected(header, storm_id, name, year, None)]

        if not matches:
            raise ValueError(
                f"No storm in {path} matches "
                f"storm_id={storm_id!r}, name={name!r}, year={year!r}.")
        if len(matches) > 1:
            found = ", ".join(f"{h['storm_id']} ({h['name']})"
                              for h in [m[0] for m in matches[:5]])
            raise ValueError(
                f"{len(matches)} storms in {path} match "
                f"storm_id={storm_id!r}, name={name!r}, year={year!r} "
                f"(e.g. {found}). Narrow the selection, or use iter_hurdat() to "
                f"read them all.")

        header, lines = matches[0]
        self = _track_from_hurdat_block(cls, header, lines, verbose=verbose)

        # HURDAT2 never supplies RMW or ROCI, so warn once for the storm rather
        # than once per track point.
        warnings.warn(missing_data_warning_str)

        self.file_paths.append(path)
        self.file_format = "hurdat"

        return self

    @classmethod
    def read_ibtracs(cls, path, sid=None, storm_name=None, year=None,
                     start_date=None, agency_pref=None):
        r"""Read in IBTrACS formatted storm file

        This reads in the netcdf-formatted IBTrACS v4 data. You must either pass
        the *sid* of the storm (a unique identifier supplied by IBTrACS) OR
        *storm_name* and *year*. The latter will not be unique for unnamed storms,
        so you may optionally pass *start_date* as well. The `wmo_\*` variable is
        used when non-missing, with missing values filled in by the corresponding
        variable of the agency specified in `wmo_agency` and/or `usa_agency`. If
        still missing, the other agencies are checked in order of *agency_pref* to
        see if any more non-missing values are available.

        :Input:
         - *path* (string) Path to the file to be read.
         - *sid* (string, optional) IBTrACS-supplied unique track identifier.
             Either *sid* OR *storm_name* and *year* must not be None.
         - *storm_name* (string, optional) name of storm of interest
             (NAME field in IBTrACS). Either *sid* OR *storm_name* and
             *year* must not be None.
         - *year* (int, optional) year of storm of interest.
             Either *sid* OR *storm_name* and *year* must not be None.
         - *start_date* (np.datetime64, optional) If storm is not
             named, will find closest unnamed storm to this start date. Only
             used for unnamed storms when specifying *storm_name* and *year*
             does not uniquely identify storm.
         - *agency_pref* (list, optional) Preference order to use if `wmo_\*` variable
             is missing and `wmo_agency` and `usa_agency` are also missing.

        :Raises:
         - *ValueError* If the method cannot find the matching storm then a
             value error is risen.
        """
        # imports that you don't need for other read functions
        try:
            import xarray as xr
        except ImportError as e:
            print("IBTrACS currently requires xarray to work.")
            raise e

        # only allow one method for specifying storms
        if (sid is not None) and ((storm_name is not None) or (year is not None)):
            raise ValueError(
                'Cannot specify both *sid* and *storm_name* or *year*.')

        if agency_pref is None:
            agency_pref = list(_IBTRACS_AGENCY_PREF)

        with xr.open_dataset(path) as ds:
            ds = _select_ibtracs_storm(ds, sid, storm_name, year, start_date)
            self = _track_from_ibtracs(cls, ds, agency_pref)

        self.file_paths.append(path)
        self.file_format = "ibtracs"

        return self


    @classmethod
    def read_jma(cls, path, verbose=False):
        r"""Read in JMA formatted storm file

        Note that only files that contain one storm are currently supported.

        For more details on the JMA format and getting data see

        http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/e_format_bst.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        :Raises:
         - *ValueError* If the method cannot find the name/year matching the
           storm or they are not provided when *single_storm == False* then a
           value error is risen.
        """
        self = cls()

        data_block = []
        with open(path, 'r') as JMA_file:
            # Extract header
            data = JMA_file.readline()
            self.ID = data[6:10]
            num_lines = int(data[12:14])
            self.name = data[30:51].strip()

            data_block = JMA_file.readlines()
        assert(num_lines == len(data_block))

        # Parse data block
        self.t = np.empty(num_lines, dtype="datetime64[s]")
        self.event = np.empty(num_lines, dtype="U2")
        self.classification = np.empty(num_lines, dtype="U4")
        self.eye_location = np.empty((num_lines, 2))
        self.max_wind_speed = np.empty(num_lines)
        self.central_pressure = np.empty(num_lines)
        self.max_wind_radius = np.empty(num_lines)
        self.storm_radius = np.empty(num_lines)
        for (i, line) in enumerate(data_block):
            if len(line) == 0:
                break
            data = [value.strip() for value in line.split()]

            # Create time from JMA yymmddhh field
            year = int(data[0][:2])
            year += 1900 if year >= 51 else 2000
            self.t[i] = np.datetime64(
                f"{year:04d}"
                f"-{data[0][2:4]}"
                f"-{data[0][4:6]}"
                f"T{data[0][6:8]}:00"
            )

            # Classification, note that this is not the category of the storm
            self.classification[i] = int(data[1])

            # Parse eye location - Always N latitude and E longitude
            self.eye_location[i, 0] = float(data[4]) / 10.0
            self.eye_location[i, 1] = float(data[3]) / 10.0

            # Intensity information - current the radii are not directly given
            # Available data includes max/min of radius of winds of 50 and
            # 30 kts instead
            self.central_pressure[i] = units.convert(
                float(data[5]), 'hPa', 'Pa')
            self.max_wind_speed[i] = units.convert(
                float(data[6]), 'knots', 'm/s')
            # Note: JMA's own missing-value convention for the numeric fields is
            # not documented in the sample data available here, so pressure and
            # wind are deliberately left unnormalized; only the radii, which the
            # format simply does not carry, are marked missing.
            self.max_wind_radius[i] = np.nan
            self.storm_radius[i] = np.nan

        # The format does not carry RMW or ROCI; warn once per file.
        warnings.warn(missing_data_warning_str)

        self.file_paths.append(path)
        self.file_format = "jma"

        return self

    @classmethod
    def read_imd(cls, path, verbose=False):
        r"""Extract relevant hurricane data from IMD file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not IMD.
        """
        raise NotImplementedError(("Reading in IMD files is not ",
                                   "implemented yet but is planned for a ",
                                   "future release."))

    @classmethod
    def read_tcvitals(cls, path, verbose=False):
        r"""Extract relevant hurricane data from TCVITALS file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        """
        self = cls()

        # read in TCVitals_file
        data_block = []
        with open(path, 'r') as TCVitals_file:
            data = TCVitals_file.readlines()
            for line in data:
                line = line.split()
                line = [value.strip() for value in line]
                data_block.append(line)
        num_lines = len(data_block)

        # Parse data block - convert to correct units
        # Conversions:
        #  max_wind_radius  - convert from km to m - 1000.0
        #  Central_pressure - convert from mbar to Pa - 100.0
        #  Radius of last isobar contour - convert from km to m - 1000.0
        self.t = np.empty(num_lines, dtype="datetime64[s]")
        self.classification = np.empty(num_lines, dtype="U4")
        self.eye_location = np.empty((num_lines, 2))
        self.max_wind_speed = np.empty(num_lines)
        self.central_pressure = np.empty(num_lines)
        self.max_wind_radius = np.empty(num_lines)
        self.storm_radius = np.empty(num_lines)

        for (i, data) in enumerate(data_block):
            # End at an empty lines - skips lines at the bottom of a file
            if len(data) == 0:
                break

            # Grab data regarding basin and cyclone number if we are starting
            if i == 0:
                self.basin = TCVitals_Basins[data[1][2:]]
                self.ID = int(data[1][:2])

            # Create time from YYYYMMDD and HHMM fields
            self.t[i] = np.datetime64(
                f"{data[3][:4]}"
                f"-{data[3][4:6]}"
                f"-{data[3][6:8]}"
                f"T{data[4][:2]}:{data[4][2:]}"
            )

            # Parse eye location - longitude/latitude order
            if data[5][-1] == 'N':
                self.eye_location[i, 1] = float(data[5][0:-1])/10.0
            else:
                self.eye_location[i, 1] = -float(data[5][0:-1])/10.0
            if data[6][-1] == "E":
                self.eye_location[i, 0] = float(data[6][0:-1])/10.0
            else:
                self.eye_location[i, 0] = -float(data[6][0:-1])/10.0

            # Intensity Information.  tcvitals marks unknown fields -99/-999;
            # normalize before unit conversion so a sentinel is never scaled
            # into a physical-looking value.
            self.max_wind_speed[i] = _sentinel_to_nan(data[12],
                                                      TCVITALS_SENTINELS)
            self.central_pressure[i] = units.convert(
                _sentinel_to_nan(data[9], TCVITALS_SENTINELS), 'mbar', 'Pa')
            self.max_wind_radius[i] = units.convert(
                _sentinel_to_nan(data[13], TCVITALS_SENTINELS), 'km', 'm')
            self.storm_radius[i] = units.convert(
                _sentinel_to_nan(data[11], TCVITALS_SENTINELS), 'km', 'm')

        self.file_paths.append(path)
        self.file_format = "tcvitals"

        return self

    # ================
    #  Track Plotting
    # ================
    def plot(self, ax, *args, t_range=None, categorization=None,
                       cat_colors={}, plot_swath=False, radius=None,
                       coordinate_system=2, fill_alpha=0.25, fill_color='red',
                       **kwargs):
        """Plot this storm's track in the given axes object

        :Input:
         - *ax* (matplotlib.pyplot.axes) Axes to plot into.
         - *t_range* (list) Time range to plot the track for.  If None then use
            entire range.  Default is None.
         - *categorization* (str) Type of categorization to be used.  This is
            used to map to the keys in the cat_colors dictionary.  Default is
            None and will cause no categorization to occur.
         - *cat_colors* (dict) Color mapping between numeric categorization and
            colors to be plotted for the track.
         - *plot_swath* (bool) Plot a swath around the track using one of the
            methods determined by what radius information is provided.  Default
            is False.
         - *radius* (None or float or numpy.ndarray)
         - *coordinate_system* (int)
         - *fill_alpha* (float)
         - *fill_color* (color)
         - *kwargs* All additional keyword arguments are passed to the plotting
            command for the track.
        """

        import matplotlib.pyplot as plt

        # Extract information for plotting the track/swath
        t = self.t
        x = self.eye_location[:, 0]
        y = self.eye_location[:, 1]
        if t_range is not None:
            t = np.ma.masked_outside(t, t_range[0], t_range[1])
            x = np.ma.array(x, mask=t.mask).compressed()
            y = np.ma.array(y, mask=t.mask).compressed()
            t = t.compressed()

        # Plot track
        if categorization is None:
            # Plot the track as a simple line with the given style
            ax.plot(x, y, *args, **kwargs)
        else:
            if self.max_wind_speed is None:
                raise ValueError("Maximum wind speed not available so "
                                 "plotting catgories is not available.")

            # Plot the track using the colors provided in the dictionary
            cat_color_defaults = {5: 'red', 4: 'yellow', 3: 'orange',
                                  2: 'green', 1: 'blue', 0: 'gray',
                                  -1: 'lightgray'}
            colors = [cat_colors.get(category, cat_color_defaults[category])
                      for category in self.category(categorization=categorization)]
            # Remove color from kwargs if they were given
            kwargs.pop('color', None)
            for i in range(t.shape[0] - 1):
                ax.plot(x[i:i+2], y[i:i+2], color=colors[i], **kwargs)

        # Plot swath
        if plot_swath:
            if (isinstance(radius, float) or isinstance(radius, np.ndarray)
                    or radius is None):

                if radius is None:
                    # Default behavior
                    if self.storm_radius is None:
                        raise ValueError("Cannot use storm radius for plotting "
                                         "the swath as the data is not available.")
                    else:
                        if coordinate_system == 1:
                            _radius = self.storm_radius
                        elif coordinate_system == 2:
                            _radius = units.convert(self.storm_radius,
                                                    'm', 'lat-long')
                        else:
                            raise ValueError(f"Unknown coordinate system "
                                             f"{coordinate_system} provided.")

                elif isinstance(radius, float):
                    # Only one value for the radius was given, replicate
                    _radius = np.ones(self.t.shape) * radius
                elif isinstance(radius, np.ndarray):
                    # The array passed is the array to use
                    _radius = radius
                else:
                    raise ValueError("Invalid input argument for radius.  Should "
                                     "be a float or None")

                # A radius that is missing or non-positive cannot be drawn --
                # matplotlib would either raise or render a nonsense patch.  The
                # readers legitimately produce NaN here (HURDAT2 carries no
                # outer radius at all), so mask those segments out rather than
                # refusing to plot the track.
                _radius = np.asarray(_radius, dtype=float)
                drawable = np.isfinite(_radius) & (_radius > 0.0)
                if not drawable.any():
                    warnings.warn(
                        "No finite positive storm radius is available, so the "
                        "swath is not drawn; pass `radius=` to plot one.")

                # Draw first and last points
                if drawable[0]:
                    ax.add_patch(plt.Circle(
                        (x[0], y[0]), _radius[0], color=fill_color))
                if t.shape[0] > 1 and drawable[-1]:
                    ax.add_patch(plt.Circle((x[-1], y[-1]), _radius[-1],
                                            color=fill_color))

                # Draw path around inner points
                if t.shape[0] > 2:
                    for i in range(t.shape[0] - 1):
                        if not (drawable[i] and drawable[i + 1]):
                            continue
                        p = np.array([(x[i], y[i]), (x[i + 1], y[i + 1])])
                        v = p[1] - p[0]
                        if abs(v[1]) > 1e-16:
                            n = np.array([1, -v[0] / v[1]], dtype=float)
                        elif abs(v[0]) > 1e-16:
                            n = np.array([-v[1] / v[0], 1], dtype=float)
                        else:
                            continue
                            # raise Exception("Zero-vector given")
                        n /= np.linalg.norm(n)
                        n *= _radius[i]

                        ax.fill((p[0, 0] + n[0], p[0, 0] - n[0],
                                 p[1, 0] - n[0],
                                 p[1, 0] + n[0]),
                                (p[0, 1] + n[1], p[0, 1] - n[1],
                                 p[1, 1] - n[1],
                                 p[1, 1] + n[1]),
                                facecolor=fill_color, alpha=fill_alpha)
                        ax.add_patch(plt.Circle((p[1][0], p[1, 1]), _radius[i],
                                                color=fill_color, alpha=fill_alpha))

    # =========================================================================
    # Other Useful Routines
    def category(self, categorization="NHC", cat_names=False):
        r"""Categorizes storm based on relevant storm data

        :Input:
         - *categorization* (string) Type of categorization to use.  Defaults
           to the National Hurricane Center "NHC".
         - *cat_names* (bool) If True returns the category name rather than a
           number.  Default to *False*.

        :Output:
         - (ndarray) Integer array of categories at each time point of the
           storm.
         - (list) Similar to the above but the name of the category as a
           *string*.  This is only returned if *car_names = True*.

        """

        # TODO:  Need to standardize on 1-minute (almost never available) or
        # 10-minute (widely available) - see
        # https://en.wikipedia.org/wiki/Tropical_cyclone#Major_basins_and_related_warning_centers

        if categorization.upper() == "BEAUFORT":
            # Beaufort scale below uses knots
            speeds = units.convert(self.max_wind_speed, "m/s", "knots")
            category = (np.zeros(speeds.shape) +
                        (speeds >= 1) * (speeds < 4) * 1 +
                        (speeds >= 4) * (speeds < 7) * 2 +
                        (speeds >= 7) * (speeds < 11) * 3 +
                        (speeds >= 11) * (speeds < 17) * 4 +
                        (speeds >= 17) * (speeds < 22) * 5 +
                        (speeds >= 22) * (speeds < 28) * 6 +
                        (speeds >= 28) * (speeds < 34) * 7 +
                        (speeds >= 34) * (speeds < 41) * 8 +
                        (speeds >= 41) * (speeds < 48) * 9 +
                        (speeds >= 48) * (speeds < 56) * 10 +
                        (speeds >= 56) * (speeds < 64) * 11 +
                        (speeds >= 64) * 12)
            cat_map = {0: "Calm",
                       1: "Light air",
                       2: "Light breeze",
                       3: "Gentle breeze",
                       4: "Moderate breeze",
                       5: "Fresh breeze",
                       6: "Strong breeze",
                       7: "High wind",
                       8: "Gale",
                       9: "Strong gale",
                       10: "Whole gale",
                       11: "Violent storm",
                       12: "Hurricane"}

        elif categorization.upper() == "NHC":
            # NHC uses knots
            speeds = units.convert(self.max_wind_speed, "m/s", "knots")
            category = (np.zeros(speeds.shape) +
                        (speeds < 34) * -1 +
                        (speeds >= 64) * (speeds < 83) * 1 +
                        (speeds >= 83) * (speeds < 96) * 2 +
                        (speeds >= 96) * (speeds < 113) * 3 +
                        (speeds >= 113) * (speeds < 137) * 4 +
                        (speeds >= 137) * 5)
            cat_map = {-1: "Tropical Depression",
                       0: "Tropical Storm",
                       1: "Category 1 Hurricane",
                       2: "Category 2 Hurricane",
                       3: "Category 3 Hurricane",
                       4: "Category 4 Hurricane",
                       5: "Category 5 Hurricane"}

        elif categorization.upper() == "JTWC":
            raise NotImplementedError("JTWC categorization not implemented.")
        elif categorization.upper() == "JMA":
            raise NotImplementedError("JMA categorization not implemented.")
        elif categorization.upper() == "IMD":
            raise NotImplementedError("IMD categorization not implemented.")
        elif categorization.upper() == "MF":
            raise NotImplementedError("MF categorization not implemented.")
        elif categorization.upper() == "BOM":
            raise NotImplementedError("BOM categorization not implemented.")
        else:
            raise ValueError("Categorization %s not available."
                             % categorization)

        if cat_names:
            category_name = []
            for (i, cat) in enumerate(category):
                category_name.append(cat_map[cat])

            return category, category_name
        else:
            return category


def _select_ibtracs_storm(ds, sid=None, storm_name=None, year=None,
                          start_date=None):
    r"""Narrow an IBTrACS dataset to the one storm requested.

    Returns the dataset with the ``storm`` dimension removed.
    """
    # match on sid
    if sid is not None:
        match = ds.sid == sid.encode()
    # or match on storm_name and year
    else:
        storm_name = storm_name.upper()
        # in case storm is unnamed
        if storm_name.upper() in ['UNNAMED', 'NO-NAME']:
            storm_name = 'NOT_NAMED'
        storm_match = (ds.name == storm_name.encode())
        year_match = (ds.time.dt.year == year).any(dim='date_time')
        match = storm_match & year_match
    # Squeeze only the storm dimension: a bare squeeze() would also
    # collapse date_time for a storm with a single observation.
    ds = ds.sel(storm=match)
    if ds.sizes.get('storm') == 1:
        ds = ds.squeeze(dim='storm')

    # occurs if we have 0 or >1 matching storms
    if 'storm' in ds.sizes:
        if ds.storm.shape[0] == 0:
            raise ValueError('Storm/year not found in provided file')
        else:
            # see if a date was provided for multiple unnamed storms
            assert start_date is not None, ValueError(
                'Multiple storms identified and no start_date specified.')

            start_times = ds.time.isel(date_time=0)
            start_date = np.datetime64(start_date)

            # find storm with start date closest to provided
            storm_ix = abs(start_times - start_date).argmin()
            ds = ds.isel(storm=storm_ix)
            assert 'storm' not in ds.sizes

    return ds


def _track_from_ibtracs(cls, ds, agency_pref):
    r"""Build a track from a single-storm IBTrACS dataset (no ``storm`` dim)."""
    import xarray as xr

    self = cls()

    # cut down dataset to only non-null times
    valid_t = ds.time.notnull()
    if valid_t.sum() == 0:
        raise ValueError('No valid wind speeds found for this storm.')
    ds = ds.sel(date_time=valid_t)

    # list of the agencies that correspond to 'usa_*' variables
    usa_agencies = [b'atcf', b'hurdat_atl', b'hurdat_epa', b'jtwc_ep',
                    b'nhc_working_bt', b'tcvightals', b'tcvitals']

    # Create mapping from wmo_ or usa_agency
    # to the appropriate variable
    agency_map = {b'': agency_pref.index('wmo')}
    # account for multiple usa agencies
    for a in usa_agencies:
        agency_map[a] = agency_pref.index('usa')
    # map all other agencies to themselves
    for i in [a for a in agency_pref if a not in ['wmo', 'usa']]:
        agency_map[i.encode('utf-8')] = agency_pref.index(i)

    # fill in usa as provider if usa_agency is
    # non-null when wmo_agency is null
    provider = ds.wmo_agency.where(ds.wmo_agency != b'', ds.usa_agency)

    # get index into from agency that is wmo_provider
    def map_val_to_ix(a):
        def func(x): return agency_map[x]
        return xr.apply_ufunc(func, a, vectorize=True)
    pref_agency_ix = map_val_to_ix(provider)

    # GET MAX WIND SPEED and PRES
    pref_vals = {}
    for v in ['wind', 'pres']:
        all_vals = ds[['{}_{}'.format(i, v) for i in agency_pref]].to_array(
            dim='agency')

        # get wmo value
        val_pref = ds['wmo_'+v]

        # fill this value in as a second-best
        pref_2 = all_vals.isel(agency=pref_agency_ix)
        val_pref = val_pref.fillna(pref_2)

        # now use the agency_pref order to fill in
        # any remaining values as third best
        best_ix = all_vals.notnull().argmax(dim='agency')
        pref_3 = all_vals.isel(agency=best_ix)
        val_pref = val_pref.fillna(pref_3)

        # add to dict
        pref_vals[v] = val_pref

    # THESE CANNOT BE MISSING SO DROP
    # IF EITHER MISSING
    valid = pref_vals['wind'].notnull() & pref_vals['pres'].notnull()
    if not valid.any():
        raise NoDataError(missing_necessary_data_warning_str)
    ds = ds.sel(date_time=valid)
    for i in ['wind', 'pres']:
        pref_vals[i] = pref_vals[i].sel(date_time=valid)

    # GET RMW and ROCI
    # (these can be missing)
    for r in ['rmw', 'roci']:
        order = ['{}_{}'.format(i, r) for i in agency_pref if
                 '{}_{}'.format(i, r) in ds.data_vars.keys()]
        vals = ds[order].to_array(dim='agency')
        best_ix = vals.notnull().argmax(dim='agency')
        val_pref = vals.isel(agency=best_ix)
        pref_vals[r] = val_pref

    # GET QUADRANT WIND RADII (also optional)
    # IBTrACS names these `<agency>_r34` / `_r50` / `_r64` along a `quadrant`
    # dimension ordered NE, SE, SW, NW -- the same order ATCF and HURDAT2 use.
    # Missing stays NaN.  Trimmed subsets of the archive omit these variables
    # entirely, so the array is left all-NaN rather than requiring them.
    quadrant_radii = np.full((ds.sizes['date_time'],
                              len(WIND_RADII_THRESHOLDS_KT), 4), np.nan)
    for j, threshold in enumerate((34, 50, 64)):
        names = ['{}_r{}'.format(agency, threshold) for agency in agency_pref
                 if '{}_r{}'.format(agency, threshold) in ds.data_vars.keys()]
        if not names or 'quadrant' not in ds.sizes:
            continue
        vals = ds[names].to_array(dim='agency')
        # Prefer the first agency reporting anything for this threshold.
        best_ix = vals.notnull().any(dim='quadrant').argmax(dim='agency')
        selected = vals.isel(agency=best_ix)
        quadrant_radii[:, j, :] = units.convert(
            selected.transpose('date_time', 'quadrant').values, 'nmi', 'm')

    # CONVERT TO GEOCLAW FORMAT

    # assign basin to be the basin where track originates
    # in case track moves across basins
    self.basin = ds.basin.values[0].astype(str)
    self.basin_code = str(self.basin)
    self.name = ds.name.astype(str).item()
    self.ID = ds.sid.astype(str).item()

    # Times as a plain datetime64 ndarray, matching every other reader.
    # Leaving this as a DataArray leaks xarray into consumers -- indexing
    # it yields 0-d DataArrays, which is what makes
    # fill_rad_w_other_source's interp fail.  Truncating to seconds also
    # drops the sub-second float roundoff IBTrACS carries in its times
    # (e.g. ...T06:00:00.000039936).
    self.t = ds.time.values.astype('datetime64[s]')

    # events
    self.event = ds.usa_record.values.astype(str)

    # time offset
    if (self.event == 'L').any():
        # if landfall, use last landfall
        self.time_offset = self.t[self.event == 'L'][-1]
    else:
        # if no landfall, use last time of storm
        self.time_offset = self.t[-1]

    # Classification, note that this is not the category of the storm.
    # Decode from the netCDF's bytes (b'TD') to str so it compares
    # against ordinary string literals.
    self.classification = ds.usa_status.values.astype(str)
    self.eye_location = np.array([ds.lon, ds.lat]).T

    # quadrant_radii was built from `ds` after it was narrowed to the times with
    # both a wind and a pressure, so it is already the right length.
    self.wind_radii = quadrant_radii
    self.wind_radii_thresholds = units.convert(
        np.array(WIND_RADII_THRESHOLDS_KT), 'knots', 'm/s')

    # Intensity information - for now, including only common, basic intensity
    # info.
    # TODO: add more detailed info for storms that have it
    #
    # xarray already carries missing values as NaN, which is the
    # in-memory contract, so these are converted straight through; the
    # former .where(..., -1) wrappers replaced NaN with a sentinel that
    # write_geoclaw's NaN-based fill/skip logic could not see.
    self.max_wind_speed = units.convert(
        pref_vals['wind'], 'knots', 'm/s').values
    self.central_pressure = units.convert(
        pref_vals['pres'], 'mbar', 'Pa').values
    self.max_wind_radius = units.convert(
        pref_vals['rmw'], 'nmi', 'm').values
    self.storm_radius = units.convert(
        pref_vals['roci'], 'nmi', 'm').values

    # warn if you have any missing vals for RMW or ROCI.  Note this is
    # deliberately `.any()`: the old `.max() == -1` test only fired when
    # *every* value was missing, so a partially-observed track warned
    # nothing and then silently lost rows at write time.
    if (np.isnan(self.max_wind_radius).any()
            or np.isnan(self.storm_radius).any()):
        warnings.warn(missing_data_warning_str)

    return self


# =============================================================================
# Multi-storm archive ingestion
#
# Both HURDAT2 and IBTrACS distribute one file per basin holding every storm on
# record, so driving an ensemble means walking a file rather than opening one per
# storm.  These are module-level generators rather than a collection class: what
# a driver needs beyond iteration is an *index* (``catalog_*``, a DataFrame), not
# a container, and a rival container object would fork the one merged object
# model.  Each ``read_*`` classmethod and its ``iter_*`` counterpart share a
# single parser, so there is exactly one parsing path to verify.
# =============================================================================

# A HURDAT2 storm header is "<basin><number><year>, <name>, <count>,", where
# count is the number of best-track records that follow.  That declared count is
# free self-validation: it must match what we parse, and the counts must sum to
# the number of data lines in the file.
_HURDAT_HEADER_RE = re.compile(r"^([A-Z]{2})(\d{2})(\d{4})\s*,\s*([^,]*?)\s*,"
                               r"\s*(\d+)\s*,")


def _hurdat_blocks(path):
    r"""Yield ``(header, data_lines)`` for each storm in a HURDAT2 file.

    *header* is a dict with ``storm_id``, ``basin_code``, ``number``, ``year``,
    ``name`` and ``num_records``.  Each block is framed by its declared record
    count, and a mismatch raises with the offending line number rather than
    silently absorbing the next storm's header (which is what the single-storm
    reader used to do -- it fed ``"AL092008"`` to ``np.datetime64``).
    """
    with open(path, "r") as hurdat_file:
        lines = [line.rstrip("\n") for line in hurdat_file]

    index = 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue

        match = _HURDAT_HEADER_RE.match(lines[index])
        if match is None:
            raise ValueError(
                f"{path}:{index + 1}: expected a HURDAT2 storm header, got "
                f"{lines[index][:60]!r}")

        basin_code, number, year, name, declared = match.groups()
        header = {"storm_id": f"{basin_code}{number}{year}",
                  "basin_code": basin_code,
                  "number": int(number),
                  "year": int(year),
                  "name": name.strip(),
                  "num_records": int(declared)}

        start = index + 1
        block = [line for line in lines[start:start + header["num_records"]]
                 if line.strip()]
        if len(block) != header["num_records"]:
            raise ValueError(
                f"{path}:{start + 1}: storm {header['storm_id']} declares "
                f"{header['num_records']} records but {len(block)} were found.")

        yield header, block
        index = start + header["num_records"]


def _hurdat_selected(header, storm_id, name, year, basins):
    r"""Whether a HURDAT2 block header matches the given filters."""
    if storm_id is not None and header["storm_id"].upper() != storm_id.upper():
        return False
    if name is not None and header["name"].upper() != name.upper():
        return False
    if year is not None:
        years = (year,) if isinstance(year, int) else tuple(year)
        if header["year"] not in years:
            return False
    if basins is not None:
        codes = ((basins,) if isinstance(basins, str) else tuple(basins))
        if header["basin_code"].upper() not in {b.upper() for b in codes}:
            return False
    return True


def _track_from_hurdat_block(cls, header, block, verbose=False):
    r"""Build a :class:`StormTrack` from one HURDAT2 header plus its records."""
    self = cls()
    self.basin = ATCF_basins.get(header["basin_code"], header["basin_code"])
    self.basin_code = header["basin_code"]
    self.name = header["name"]
    # ATCF-style id ("AL092008").  This used to be the header's record count.
    self.ID = header["storm_id"]

    num_lines = len(block)
    self.t = np.empty(num_lines, dtype="datetime64[s]")
    self.event = np.empty(num_lines, dtype="U2")
    self.classification = np.empty(num_lines, dtype="U4")
    self.eye_location = np.empty((num_lines, 2))
    self.max_wind_speed = np.empty(num_lines)
    self.central_pressure = np.empty(num_lines)
    self.max_wind_radius = np.empty(num_lines)
    self.storm_radius = np.empty(num_lines)
    # (time, threshold, quadrant) quadrant wind radii, SI metres.
    self.wind_radii = np.full(
        (num_lines, len(WIND_RADII_THRESHOLDS_KT), 4), np.nan)
    self.wind_radii_thresholds = units.convert(
        np.array(WIND_RADII_THRESHOLDS_KT), 'knots', 'm/s')

    for (i, line) in enumerate(block):
        data = [value.strip() for value in line.split(",")]

        # Create time from YYYYMMDD and HHMM fields
        self.t[i] = np.datetime64(
            f"{data[0][:4]}"
            f"-{data[0][4:6]}"
            f"-{data[0][6:8]}"
            f"T{data[1][:2]}:{data[1][2:]}"
        )

        # If an event is occuring record it.  If landfall then use as an
        # offset.   Note that if there are multiple landfalls the last one
        # is used as the offset
        if len(data[2].strip()) > 0:
            self.event[i] = data[2].strip()
            if self.event[i].upper() == "L":
                self.time_offset = self.t[i]

        # Classification, note that this is not the category of the storm
        self.classification[i] = data[3]

        # Parse eye location
        if data[4][-1] == "N":
            self.eye_location[i, 1] = float(data[4][0:-1])
        else:
            self.eye_location[i, 1] = -float(data[4][0:-1])
        if data[5][-1] == "E":
            self.eye_location[i, 0] = float(data[5][0:-1])
        else:
            self.eye_location[i, 0] = -float(data[5][0:-1])

        # Intensity information - radii are not included directly in this
        # format and instead radii of winds above a threshold are included.
        # HURDAT2 marks unknown wind as -99 and unknown pressure as -999;
        # normalize before converting units, or the sentinel is scaled into a
        # physical-looking value (-999 mbar -> -99900 Pa).
        self.max_wind_speed[i] = units.convert(
            _sentinel_to_nan(data[6], HURDAT_SENTINELS), 'knots', 'm/s')
        self.central_pressure[i] = units.convert(
            _sentinel_to_nan(data[7], HURDAT_SENTINELS), 'mbar', 'Pa')

        # Quadrant wind radii: fields 9-20 are the 34/50/64-kt radii for
        # NE, SE, SW, NW.  A zero here is *real data* -- a weak system simply
        # has no 34-kt winds, so their radius is zero -- and is distinct from
        # HURDAT2's -999 "unknown", which _sentinel_to_nan turns into NaN.
        if len(data) >= 20:
            radii = _sentinel_to_nan(list(data[8:20]), HURDAT_SENTINELS)
            self.wind_radii[i, :, :] = units.convert(
                radii.reshape(3, 4), 'nmi', 'm')

        # Field 21 is the radius of maximum wind, added with the 2021 season
        # (populated for 2021 onward, -999 before).  Verified against IBTrACS
        # `usa_rmw`: 2676/2676 coincident 2021-2025 North Atlantic records agree
        # exactly when read as nautical miles.
        # Revisions before 2021 end the line at 20 fields, so treat an absent
        # field the same as -999.
        self.max_wind_radius[i] = units.convert(
            _sentinel_to_nan(data[20] if len(data) >= 21 else "",
                             HURDAT_SENTINELS), 'nmi', 'm')

        # HURDAT2 carries no outer/last-closed-isobar radius in any revision.
        self.storm_radius[i] = np.nan

    if verbose:
        print(f"  {self.ID} {self.name}: {num_lines} records "
              f"{self.t[0]} -> {self.t[-1]}")
    return self


def iter_hurdat(path, storm_ids=None, names=None, years=None, basins=None,
                verbose=False):
    r"""Iterate the storms in a HURDAT2 file as :class:`StormTrack` objects.

    One pass over the file, yielding in file order.  Every filter is optional and
    they combine; ``years`` accepts an int or any iterable of ints, so a season
    range is ``years=range(1980, 2026)``.

    :Input:
     - *path* (string) Path to a HURDAT2 file (one or many storms).
     - *storm_ids* (iterable) ATCF-style ids to keep, e.g. ``["AL092008"]``.
     - *names* (iterable) Storm names to keep (case-insensitive).
     - *years* (int or iterable) Seasons to keep.
     - *basins* (string or iterable) Two-letter basin codes to keep.
     - *verbose* (bool) Print each storm as it is read.

    :Output:
     - Generator of :class:`StormTrack`, each with ``file_format = "hurdat"``.

    Missing RMW/ROCI is warned about once for the whole iteration rather than
    once per storm; HURDAT2 never carries either.
    """
    wanted_ids = (None if storm_ids is None
                  else {str(i).upper() for i in storm_ids})
    wanted_names = (None if names is None
                    else {str(n).upper() for n in names})

    warned = False
    for header, block in _hurdat_blocks(path):
        if wanted_ids is not None and header["storm_id"].upper() not in wanted_ids:
            continue
        if wanted_names is not None and header["name"].upper() not in wanted_names:
            continue
        if not _hurdat_selected(header, None, None, years, basins):
            continue

        track = _track_from_hurdat_block(StormTrack, header, block,
                                         verbose=verbose)
        if not warned:
            warnings.warn(missing_data_warning_str)
            warned = True
        track.file_paths.append(path)
        track.file_format = "hurdat"
        yield track


def _ibtracs_seasons(ds):
    r"""Season (year) per storm, from ``season`` if present else from ``time``.

    The variable is present in the distributed basin files but absent from
    trimmed subsets, so fall back rather than requiring it.
    """
    if "season" in ds.variables:
        return ds["season"].values.astype(int)
    return ds.time.isel(date_time=0).dt.year.values.astype(int)


def _ibtracs_storm_mask(ds, basin=None, years=None, sids=None, names=None):
    r"""Boolean mask over the ``storm`` dimension for the requested filters."""
    mask = np.ones(ds.sizes["storm"], dtype=bool)

    if basin is not None:
        codes = ((basin,) if isinstance(basin, str) else tuple(basin))
        wanted = {b.upper().encode() for b in codes}
        # basin varies along a track; select on where the storm originates,
        # matching what the reader stores as `basin`.
        first = ds.basin.isel(date_time=0).values
        mask &= np.array([b.upper() in wanted for b in first])

    if years is not None:
        wanted_years = ([int(years)] if isinstance(years, (int, np.integer))
                        else [int(y) for y in years])
        mask &= np.isin(_ibtracs_seasons(ds), wanted_years)

    if sids is not None:
        wanted_sids = {str(s).upper().encode() for s in sids}
        mask &= np.array([s.upper() in wanted_sids for s in ds.sid.values])

    if names is not None:
        wanted_names = {str(n).upper().encode() for n in names}
        mask &= np.array([n.upper() in wanted_names for n in ds.name.values])

    return mask


def iter_ibtracs(path, basin=None, years=None, sids=None, names=None,
                 agency_pref=None, skip_invalid=True, verbose=False):
    r"""Iterate storms in an IBTrACS v4 netCDF file as :class:`StormTrack`s.

    Opens the file once and slices per storm, rather than re-scanning it for
    every :meth:`StormTrack.read_ibtracs` call.

    :Input:
     - *path* (string) Path to an IBTrACS v4 netCDF file.
     - *basin* (string or iterable) Two-letter basin code(s) of storm *origin*,
       e.g. ``"NA"``; matches what the reader stores as ``basin``.
     - *years* (int or iterable) Seasons to keep, e.g. ``range(1980, 2026)``.
     - *sids* (iterable) IBTrACS storm ids to keep.
     - *names* (iterable) Storm names to keep (case-insensitive).
     - *agency_pref* (list) As :meth:`StormTrack.read_ibtracs`.
     - *skip_invalid* (bool) Skip storms with no time at which both a wind and a
       pressure were reported, instead of raising.  Load-bearing over a whole
       basin: such storms exist, and one of them should not abort the sweep.
     - *verbose* (bool) Report each storm as it is read.

    :Output:
     - Generator of :class:`StormTrack`, in file order.
    """
    try:
        import xarray as xr
    except ImportError as e:
        print("IBTrACS currently requires xarray to work.")
        raise e

    if agency_pref is None:
        agency_pref = list(_IBTRACS_AGENCY_PREF)

    n_skipped = 0
    with xr.open_dataset(path) as ds:
        if "storm" not in ds.sizes:
            raise ValueError(
                f"{path} has no 'storm' dimension; it does not look like an "
                "IBTrACS basin file.")

        mask = _ibtracs_storm_mask(ds, basin=basin, years=years, sids=sids,
                                   names=names)
        indices = np.nonzero(mask)[0]
        if verbose:
            print(f"{len(indices)} of {ds.sizes['storm']} storms selected "
                  f"from {path}")

        for index in indices:
            one = ds.isel(storm=int(index))
            sid = one.sid.astype(str).item()
            try:
                with warnings.catch_warnings():
                    # Missing RMW/ROCI is the norm across a whole basin; warn
                    # once for the sweep rather than once per storm.
                    warnings.simplefilter("ignore", UserWarning)
                    track = _track_from_ibtracs(StormTrack, one, agency_pref)
            except (NoDataError, ValueError) as error:
                if not skip_invalid:
                    raise
                n_skipped += 1
                if verbose:
                    print(f"  skipped {sid}: {error}")
                continue

            track.file_paths.append(path)
            track.file_format = "ibtracs"
            if verbose:
                print(f"  {track.ID} {track.name}: {len(track.t)} records "
                      f"{track.t[0]} -> {track.t[-1]}")
            yield track

    if n_skipped:
        warnings.warn(f"Skipped {n_skipped} IBTrACS storm(s) with no time "
                      "having both a wind and a pressure observation.")


def catalog_ibtracs(path, basin=None, years=None):
    r"""Index of the storms in an IBTrACS file, as a :class:`pandas.DataFrame`.

    Columns: ``sid``, ``name``, ``season``, ``basin``, ``num_records``,
    ``t_start``, ``t_end``, ``max_wind`` (m/s, from ``wmo_wind`` where present).
    Reads summary variables only, so it is far cheaper than building every
    track -- use it to decide what to run before running it.
    """
    try:
        import xarray as xr
    except ImportError as e:
        print("IBTrACS currently requires xarray to work.")
        raise e

    with xr.open_dataset(path) as ds:
        mask = _ibtracs_storm_mask(ds, basin=basin, years=years)
        indices = np.nonzero(mask)[0]
        seasons = _ibtracs_seasons(ds)
        times = ds.time.values
        origin_basin = ds.basin.isel(date_time=0).values

        rows = []
        for index in indices:
            valid = ~np.isnat(times[index])
            if not valid.any():
                continue
            track_times = times[index][valid]
            wind = ds.wmo_wind.isel(storm=int(index)).values
            rows.append({
                "sid": ds.sid.values[index].astype(str),
                "name": ds.name.values[index].astype(str),
                "season": int(seasons[index]),
                "basin": origin_basin[index].astype(str),
                "num_records": int(valid.sum()),
                "t_start": track_times[0].astype("datetime64[s]"),
                "t_end": track_times[-1].astype("datetime64[s]"),
                "max_wind": (float(units.convert(np.nanmax(wind), 'knots', 'm/s'))
                             if np.isfinite(wind).any() else np.nan),
            })

    return pd.DataFrame(rows, columns=["sid", "name", "season", "basin",
                                       "num_records", "t_start", "t_end",
                                       "max_wind"])



def catalog_hurdat(path, years=None, basins=None):
    r"""Index of the storms in a HURDAT2 file, as a :class:`pandas.DataFrame`.

    Columns: ``storm_id``, ``name``, ``year``, ``basin``, ``num_records``,
    ``t_start``, ``t_end``.  Reads only each storm's first and last record, so it
    is much cheaper than building every track -- use it to decide what to run
    before running it.
    """
    rows = []
    for header, block in _hurdat_blocks(path):
        if not _hurdat_selected(header, None, None, years, basins):
            continue

        def _time(line):
            fields = [value.strip() for value in line.split(",")]
            return np.datetime64(f"{fields[0][:4]}-{fields[0][4:6]}"
                                 f"-{fields[0][6:8]}T{fields[1][:2]}"
                                 f":{fields[1][2:]}")

        rows.append({"storm_id": header["storm_id"],
                     "name": header["name"],
                     "year": header["year"],
                     "basin": ATCF_basins.get(header["basin_code"],
                                              header["basin_code"]),
                     "num_records": header["num_records"],
                     "t_start": _time(block[0]),
                     "t_end": _time(block[-1])})
    return pd.DataFrame(rows, columns=["storm_id", "name", "year", "basin",
                                       "num_records", "t_start", "t_end"])


# =============================================================================
# Radius fill functions
def fill_rad_w_other_source(t, storm_targ, storm_fill, var, interp_kwargs={}):
    r"""Fill in storm radius variable (*max_wind_radius* or \
    *storm_radius*) with values from another source. i.e.
    if you have missing radii in IBTrACS, you can fill with ATCF.
    This function will assume *storm_fill* has more non-missing
    values than *storm_targ* for this particular radius variable.
    Thus, it first attempts to interpolate the variable in *storm_fill*
    to the desired timestep. If that is missing, it tries to interpolate
    the non-missing values of the variable in *storm_targ*. If that
    also fails, it returns ``np.nan``. The proper usage of this
    function is to wrap it such that you can pass a function
    with (*t*, *storm*) arguments to *max_wind_radius_fill* or
    *storm_radius_fill* when calling *write_geoclaw*.

    :Input:
     - *t* (np.datetime64) the time corresponding to
       a missing value of *max_wind_radius* or *storm_radius*
     - *storm_targ* (:py:class:`clawpack.geoclaw.storm.Storm`) storm
       that has missing values you want to fill
     - *storm_fill* (:py:class:`clawpack.geoclaw.storm.Storm`) storm
       that has non-missing values you want to use to fill *storm_targ*
     - *var* (str) Either 'max_wind_radius' or 'storm_radius'
     - *interp_kwargs* (dict) Additional keywords passed to scipy's
       interpolator.

    :Returns:
     - (float) value to use to fill this time point in *storm_targ*.
       ``np.nan`` if still missing after using *storm_fill* to fill, so a
       failed fill is skipped by ``write_geoclaw`` rather than written out.

    :Examples:

    .. code-block:: python

        >>> storm_ibtracs = Storm(file_format='IBTrACS', path='path_to_ibtracs.nc',
        ...     sid='2018300N26315')

        >>> storm_atcf = Storm(file_format='ATCF', path='path_to_atcf.dat')

        >>> def fill_mwr(t, storm):
        ...     return fill_rad_w_other_source(t, storm, storm_atcf, 'max_wind_radius')

        >>> storm_ibtracs.write(file_format = 'geoclaw',
        ...     path = 'out_path.storm',
        ...     max_wind_radius_fill = fill_mwr)
    """

    try:
        import xarray as xr
    except ImportError as e:
        print("fill_rad_w_other_source currently requires xarray to work.")
        raise e

    # Accept a 0-d DataArray for *t* as well as a scalar.  Callers commonly pass
    # ``storm.t[n]``, and a reader that stores its times as a DataArray hands back
    # a 0-d DataArray, which interp() below cannot convert to a datetime.
    if hasattr(t, 'values'):
        t = t.values
    if isinstance(t, np.ndarray) and t.ndim == 0:
        t = t.item()

    fill_da = xr.DataArray(getattr(storm_fill, var),
                           coords={'t': np.asarray(getattr(storm_fill, 't'))},
                           dims=('t',))

    # convert -1 to nan
    fill_da = fill_da.where(fill_da > 0, np.nan)

    # if not all missing, try using storm_fill to fill
    if fill_da.notnull().any():

        # remove duplicates
        fill_da = fill_da.groupby('t').first()

        # remove NaNs
        fill_da = fill_da.dropna('t')

        # interpolate to point
        fill_interp = fill_da.interp(t=[t], kwargs=interp_kwargs).item()

        # try replacing with storm_fill
        # (assuming atcf has more data points than ibtracs)
        if not np.isnan(fill_interp):
            return fill_interp

    # next, try just interpolating other ibtracs values
    targ_da = xr.DataArray(getattr(storm_targ, var),
                           coords={'t': np.asarray(getattr(storm_targ, 't'))},
                           dims=('t',))
    targ_da = targ_da.where(targ_da > 0, np.nan)
    if targ_da.notnull().any():
        targ_da = targ_da.groupby('t').first()
        targ_da = targ_da.dropna('t')
        targ_interp = targ_da.interp(t=[t], kwargs=interp_kwargs).item()
        if not np.isnan(targ_interp):
            return targ_interp

    # If nothing worked, report missing.  NaN rather than -1 so write_geoclaw
    # re-checks the fill's return and skips the row instead of emitting a
    # negative radius that GeoClaw cannot run.  The `> 0` masks above are kept:
    # for a radius that is the physically right predicate, and NaN > 0 is False.
    return np.nan
