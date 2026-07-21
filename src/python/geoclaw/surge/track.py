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
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.max_wind_speed = None
        self.max_wind_radius = None
        self.central_pressure = None
        self.storm_radius = None
        self.classification = None
        self.basin = None
        self.wind_speeds = None

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
        self.basin = ATCF_basins[df["BASIN"][0]]
        self.ID = df["CY"][0]

        # Keep around the name as an array
        self.name = df["STORMNAME"].to_numpy()

        # Take forecast period TAU into consideration
        df['DATE'] = df["YYYYMMDDHH"] + df["TAU"]
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

        self.file_paths.append(path)
        self.file_format = "atcf"

        return self

    @classmethod
    def read_hurdat(cls, path, verbose=False):
        r"""Read in HURDAT formatted storm file

        This is the current version of HURDAT data available (HURDAT 2).  Note
        that this assumes there is only one storm in the file (includes the
        header information though).  Future features will be added that will allow for
        a file to be read with multiple storms defined.

        For more details on the HURDAT format and getting data see

        http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        :Raises:
         - *ValueError* If the method cannot find the name/year matching the
           storm or they are not provided when *single_storm == False* then a
           value error is risen.
        """
        self = cls()

        with open(path, 'r') as hurdat_file:
            # Extract header
            data = [value.strip() for value in
                    hurdat_file.readline().split(',')]
            self.basin = data[0][:2]
            self.name = data[1]
            self.ID = data[2]

            # Store rest of data
            data_block = hurdat_file.readlines()

        num_lines = len(data_block)

        # Parse data block
        self.t = np.empty(num_lines, dtype="datetime64[s]")
        self.event = np.empty(num_lines, dtype=str)
        self.classification = np.empty(num_lines, dtype=str)
        self.eye_location = np.empty((num_lines, 2))
        self.max_wind_speed = np.empty(num_lines)
        self.central_pressure = np.empty(num_lines)
        self.max_wind_radius = np.empty(num_lines)
        self.storm_radius = np.empty(num_lines)

        for (i, line) in enumerate(data_block):
            if len(line) == 0:
                break
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
            # format and instead radii of winds above a threshold are included
            self.max_wind_speed[i] = units.convert(
                float(data[6]), 'knots', 'm/s')
            self.central_pressure[i] = units.convert(
                float(data[7]), 'mbar', 'Pa')
            warnings.warn(missing_data_warning_str)
            self.max_wind_radius[i] = -1
            self.storm_radius[i] = -1

        self.file_paths.append(path)
        self.file_format = "hurdat"

        return self

    @classmethod
    def read_ibtracs(cls, path, sid=None, storm_name=None, year=None, start_date=None,
                     agency_pref=['wmo',
                                  'usa',
                                  'tokyo',
                                  'newdelhi',
                                  'reunion',
                                  'bom',
                                  'nadi',
                                  'wellington',
                                  'cma',
                                  'hko',
                                  'ds824',
                                  'td9636',
                                  'td9635',
                                  'neumann',
                                  'mlc']):
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
        self = cls()

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

        with xr.open_dataset(path) as ds:

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
            ds = ds.sel(storm=match).squeeze()

            # occurs if we have 0 or >1 matching storms
            if 'storm' in ds.dims.keys():
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
                    ds = ds.isel(storm=storm_ix).squeeze()
                    assert 'storm' not in ds.dims.keys()

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

            # CONVERT TO GEOCLAW FORMAT

            # assign basin to be the basin where track originates
            # in case track moves across basins
            self.basin = ds.basin.values[0].astype(str)
            self.name = ds.name.astype(str).item()
            self.ID = ds.sid.astype(str).item()

            # convert datetime64 to datetime.datetime
            self.t = ds.time

            # events
            self.event = ds.usa_record.values.astype(str)

            # time offset
            if (self.event == 'L').any():
                # if landfall, use last landfall
                self.time_offset = np.array(self.t)[self.event == 'L'][-1]
            else:
                # if no landfall, use last time of storm
                self.time_offset = self.t[-1]

            # Classification, note that this is not the category of the storm
            self.classification = ds.usa_status.values
            self.eye_location = np.array([ds.lon, ds.lat]).T

            # Intensity information - for now, including only common, basic intensity
            # info.
            # TODO: add more detailed info for storms that have it
            self.max_wind_speed = units.convert(
                pref_vals['wind'], 'knots', 'm/s').where(pref_vals['wind'].notnull(), -1).values
            self.central_pressure = units.convert(pref_vals['pres'], 'mbar', 'Pa').where(
                pref_vals['pres'].notnull(), -1).values
            self.max_wind_radius = units.convert(pref_vals['rmw'], 'nmi', 'm').where(
                pref_vals['rmw'].notnull(), -1).values
            self.storm_radius = units.convert(pref_vals['roci'], 'nmi', 'm').where(
                pref_vals['roci'].notnull(), -1).values

            # warn if you have missing vals for RMW or ROCI
            if (self.max_wind_radius.max()) == -1 or (self.storm_radius.max() == -1):
                warnings.warn(missing_data_warning_str)

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
        self.event = np.empty(num_lines, dtype=str)
        self.classification = np.empty(num_lines, dtype=str)
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
            warnings.warn(missing_data_warning_str)
            self.max_wind_radius[i] = -1
            self.storm_radius[i] = -1

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
        self.classification = np.empty(num_lines, dtype=str)
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

            # Intensity Information
            self.max_wind_speed[i] = float(data[12])
            self.central_pressure[i] = units.convert(
                float(data[9]), 'mbar', 'Pa')
            self.max_wind_radius[i] = units.convert(float(data[13]), 'km', 'm')
            self.storm_radius[i] = units.convert(float(data[11]), 'km', 'm')

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

                # Draw first and last points
                ax.add_patch(plt.Circle(
                    (x[0], y[0]), _radius[0], color=fill_color))
                if t.shape[0] > 1:
                    ax.add_patch(plt.Circle((x[-1], y[-1]), _radius[-1],
                                            color=fill_color))

                # Draw path around inner points
                if t.shape[0] > 2:
                    for i in range(t.shape[0] - 1):
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
    also fails, it simply returns -1. The proper usage of this
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
    - (float) value to use to fill this time point in *storm_targ*. -1
        if still missing after using *storm_fill* to fill.

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

    # if nothing worked, return the missing value (-1)
    return -1
