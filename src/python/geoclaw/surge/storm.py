#!/usr/bin/env python

r"""
Module defines a class and routines for managing storm best-track type input and
testing reconstructed wind and pressure fields.  Additionally some support for
ensembles of storms from various providers is also included.

The primary class of interest in the module is the `Storm` class that
facilitates dealing with various best-track formats found around the world and
the expected GeoClaw storm format that is read into the FORTRAN code.  The basic
workflow in a `setrun.py` file would do the following:

1. Create a `Storm` object by reading in from a file::

    storm = clawpack.geoclaw.surge.storm.Storm("my_storm.txt", file_format='ATCF')

2. Write out the storm object created into the GeoClaw format::

    storm.write("my_geoclaw_storm.txt", file_format="geoclaw")

3. Specify the path to the GeoClaw formatted storm file, in this case
   "my_geoclaw_storm.txt".

:Formats Supported:
    - GeoClaw (fully)
    - ATCF (reading only)
    - HURDAT (reading only)
    - IBTrACS (reading only)
    - JMA (reading only)
    - IMD (planned)
    - tcvitals (reading only)
    - Data
"""

import sys
# import os
import argparse
import datetime
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import clawpack.geoclaw.units as units
import clawpack.geoclaw.util as util
import clawpack.clawutil.data as clawdata


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


# =============================================================================
#  Basic storm class
class Storm(object):
    r"""
    Storm data object

    This object contains a time series of time data that describe a particular
    storm.  This includes the attributes below and the ability to read from
    multiple sources for data such as the U.S. National Hurricane Center (NHC),
    the Japanese Meterological Agency (JMA), and the Indian Meteorlogical
    Department (IMD).  This class can then write out in any of these formats,
    construct the wind and pressure fields using a supported parameterized
    model, or output the GeoClaw supported storm format used for running storm
    surge simulations.

    *TODO:*  Add description of unit handling

    :Attributes:
     - *t* (list(float) or list(np.datetiem64)) Contains the time at which
       each entry of the other arrays are at.  These are expected to
       be *datetime64* objects. Note that when written some formats require
       a *time_offset* to be set.
     - *eye_location* (ndarray(:, :)) location of the eye of the storm. Default
       units are in signed decimal longitude and latitude.
     - *max_wind_speed* (ndarray(:)) Maximum wind speed.  Default units are
       meters/second.
     - *max_wind_radius* (ndarray(:)) Radius at which the maximum wind speed
       occurs.  Default units are meters.
     - *central_pressure* (ndarray(:)) Central pressure of storm.  Default units
       are Pascals.
     - *storm_radius* (ndarray(:)) Radius of storm, often defined as the last
       closed iso-bar of pressure.  Default units are meters.
     - *time_offset* (np.datetiem64) A date time that as an offset for the
       simulation time.  This will default to the beginning of the first of
       the year that the first time point is found in.
     - *wind_speeds* (ndarray(:, :)) Wind speeds defined in every record, such
       as 34kt, 50kt, 64kt, etc and their radii. Default units are
       meters/second and meters.

    :Initialization:
     1. Read in existing file at *path*.
     2. Construct an empty storm and supply the fields needed.  Note that these
        fields must be converted to the appropriate units.

    :Input:
     - *path* (string) Path to file to be read in if requested.
     - *file_format* (string) Format of file at path.  Default is "hurdat"
     - *kwargs* (dict) Other key-word arguments are passed to the appropriate
       read routine.
    """

    # Define supported formats and models - keys are function name related and
    # values are the proper name and a citation or URL documenting the format
    _supported_formats = {"geoclaw": ["GeoClaw", "http://www.clawpack.org/storms"],
                          "atcf": ["ATCF", "http://www.nrlmry.navy.mil/atcf_web/docs/database/new/database.html"],
                          "hurdat": ["HURDAT", "http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html"],
                          "ibtracs": ["IBTrACS", "https://www.ncdc.noaa.gov/ibtracs/index.php?name=ib-v4-access"],
                          "jma": ["JMA", "http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/e_format_bst.html"],
                          "imd": ["IMD", "http://www.rsmcnewdelhi.imd.gov.in/index.php"],
                          "tcvitals": ["TC-Vitals", "http://www.emc.ncep.noaa.gov/mmb/data_processing/tcvitals_description.htm"],
                          "netcdf": ["NetCDF", None],
                          "owi": ['OWI', "http://www.oceanweather.com"],
                          "data": ["GeoClaw Data", "http://www.clawpack.org/storms"]}


    def __init__(self, path=None, file_format="ATCF", **kwargs):
        r"""Storm Initiatlization Routine

        See :class:`Storm` for more info.
        """

        # Time offsets are usually set to landfall but could be any time point
        # and are not required
        self.time_offset = None
        # File paths of either the original file that was read in for modeled
        # storms or a list of files to be pointed to for data driven storms
        self.file_paths = []
        # Either the format that was read in or the type of data-driven storm
        self.file_format = None

        # Model parameters stored directly in the storm file
        self.t = None
        self.eye_location = None
        self.max_wind_speed = None
        self.max_wind_radius = None
        self.central_pressure = None
        self.storm_radius = None
        self.wind_speeds = None

        # Storm descriptions - not all formats provide these
        self.name = None                    # Possibly a list of a storm's names
        self.basin = None                   # Basin containing storm
        self.ID = None                      # ID code - depends on format
        self.classification = None          # Classification of storm (e.g. HU)
        self.event = None                   # Event (e.g. landfall) - HURDAT

        # Ramping information - only applies to data storms currently
        self.window_type = 0
        self.ramp_width = 1e0               # 1 degree
        self.window = None                  # If not provided for data files it
                                            # will be set to the file extents

        if path is not None:
            self.read(path, file_format=file_format, **kwargs)

    # ==========================================================================
    #  Basic object support
    def __str__(self):
        r""""""
        output = f"Name: {self.name}\n"
        if self.t is None and self.time_offset is not None:
            output += f"Time offset: {self.time_offset}\n"
        elif isinstance(self.t[0], np.datetiem64):
            output += f"Dates: {self.t[0].isoformat()}"
            output += f" - {self.t[-1].isoformat()}\n"
        else:
            output += f"Dates: {self.t[0]} - {self.t[-1]}\n"
        output += "File paths:"
        for path in self.file_paths:
            output += f"\n  {path}"

        return output

    def __repr__(self):
        return '<{}.{} "{}" at {}>'.format(
            self.__class__.__module__,
            self.__class__.__name__,
            self.__dict__.get('name', 'name not given'),
            hex(id(self)))

    # ==========================================================================
    # Read Routines
    def read(self, path=None, file_format="atcf", **kwargs):
        r"""Read in storm data from *path* with format *file_format*

        :Input:
         - *path* (string) Path to data file.
         - *file_format* (string) Format of the data file.  See list of
           supported formats for a list of valid strings.  Defaults to
           "hurdat".
         - *kwargs* (dict) Keyword dictionary for additional arguments that can
           be passed down to the appropriate read functions.  Please refer to
           the specific routine for a list of valid options.

        :Raises:
         - *ValueError* If the *file_format* requested does not match any of
           the available supported formats a *ValueError* is raised.
        """

        # If a path is not provided then we can try and find the relevant
        # database and download it
        if path is None:
            data_str = ("Currently automatic download of storm databases is ",
                        "not implemented.  Please refer to the URLs below for",
                        "references as to where you can download storm data",
                        "files:",
                        " - ATCF - http://ftp.nhc.noaa.gov/atcf/archive/",
                        " - HURDAT - http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html",
                        " - IBTrACS - https://www.ncdc.noaa.gov/ibtracs/index.php?name=ib-v4-access",
                        " - JMA - http://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/besttrack.html",
                        " - IMD - http://www.rsmcnewdelhi.imd.gov.in/index.php",
                        " - TCVITALS - http://www.emc.ncep.noaa.gov/mmb/data_processing/tcvitals_description.htm")
            raise NotImplementedError("\n".join(data_str))

        if file_format.lower() not in self._supported_formats.keys():
            raise ValueError("File format %s not available." % file_format)

        if isinstance(path, str):
            path = Path(path)
        getattr(self, 'read_%s' % file_format.lower())(path, **kwargs)

    def read_geoclaw(self, path, verbose=False):
        r"""Read in a GeoClaw formatted storm file

        GeoClaw storm files are read in by the Fortran code and are not meant
        to be human readable.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
        """

        # Attempt to get name from file if is follows the convention name.storm
        if path.suffix == ".storm":
            self.name = path.name

        # Read header
        with open(path, 'r') as data_file:
            num_casts = int(data_file.readline())
            time = data_file.readline()[:19]
            try:
                self.time_offset = np.datetime64(time)
            except ValueError:
                self.time_offset = float(time)
        # Read rest of data
        data = np.loadtxt(path, skiprows=3)
        num_forecasts = data.shape[0]
        self.eye_location = np.empty((num_forecasts, 2))
        assert(num_casts == num_forecasts)
        if isinstance(self.time_offset, np.datetime64):
            self.t = np.array([self.time_offset
                               + np.timedelta64(data[i, 0], "s")
                               for i in range(num_forecasts)])
        else:
            self.t = data[:, 0]
        self.eye_location[:, 0] = data[:, 1]
        self.eye_location[:, 1] = data[:, 2]
        self.max_wind_speed = data[:, 3]
        self.max_wind_radius = data[:, 4]
        self.central_pressure = data[:, 5]
        self.storm_radius = data[:, 6]

        self.file_paths.append(path)
        self.file_format = "geoclaw"

    def read_atcf(self, path, verbose=False):
        r"""Read in a ATCF formatted storm file

        ATCF format has storm stored individually so there is no support for
        multiple storms in a particular file.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
        """

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

    def read_hurdat(self, path, verbose=False):
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
        self.t = np.empty(num_lines, dtype=np.datetime64)
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

            # Create time
            self.t[i] = np.datetime64(f"{data[0][:4]}"      + 
                                      f"-{data[0][4:6]}"    + 
                                      f"-{data[0][6:8]}"    + 
                                      f"T{data[1][:2]}"     + 
                                      f":{data[1][2:]}")

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

    def read_ibtracs(self, path, sid=None, storm_name=None, year=None, start_date=None,
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

    def read_jma(self, path, verbose=False):
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
        self.t = np.empty(num_lines, dtype=np.datetime64)
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

            # Create time
            self.t[i] = np.datetime64(f"{data[0][:2]}"      + 
                                      f"-{data[0][2:4]}"    + 
                                      f"-{data[0][4:6]}"    + 
                                      f"T{data[0][6:]}")

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

    def read_imd(self, path, verbose=False):
        r"""Extract relevant hurricane data from IMD file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.

        Return ValueError if format incorrect or if file not IMD.
        """
        raise NotImplementedError(("Reading in IMD files is not ",
                                   "implemented yet but is planned for a ",
                                   "future release."))

    def read_tcvitals(self, path, verbose=False):
        r"""Extract relevant hurricane data from TCVITALS file
            and update storm fields with proper values.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.

        """

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
        self.t = np.empty(num_lines, dtype=np.datetime64)
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

            # Create time
            self.t[i] = np.datetime64(f"{data[0][:2]}"      + 
                                      f"-{data[0][2:4]}"    + 
                                      f"-{data[0][4:6]}"    + 
                                      f"T{data[0][6:]}")

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

    def read_data(self, path, verbose=False):
        r"""Read in a data derived storm data information file

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
        """

        with path.open() as data_file:
            data_file.readline()
            self.time_offset = np.datetime64(data_file.readline()[:19])
            self.file_format = int(
                                data_file.readline().partition("#")[0].rstrip())
            num_files = int(data_file.readline().partition("#")[0].rstrip())
            self.window_type = int(
                                data_file.readline().partition("#")[0].rstrip())
            self.ramp_width = float(
                                data_file.readline().partition("#")[0].rstrip())
            # self.window = [value for value in 
            #             data_file.readline().partition("#")[0].rstrip().split()]
            data_file.readline()
            data_file.readline()
            data_file.readline()

            # We do not keep track of any of the format specific information
            if self.file_format == 1:
                data_file.readline()
                data_file.readline()
            elif self.file_format == 2:
                data_file.readline()
                data_file.readline()
                data_file.readline()
                data_file.readline()
            else:
                raise TypeError(f"Unknown storm data file format type" +
                                f" '{self.file_format}' provided.")

            for i in range(num_files):
                self.file_paths.append(Path(data_file.readline().rstrip()))


    # =========================================================================
    # Write Routines
    def write(self, path, file_format="geoclaw", **kwargs):
        r"""Write out the storm data to *path* in format *file_format*

        :Input:
         - *path* (string) Path to data file.
         - *file_format* (string) Format of the data file.  See list of
           supported formats for a list of valid strings.  Defaults to
           "geoclaw".
         - *kwargs* (dict) Keyword dictionary for additional arguments that can
           be passed down to the appropriate write functions.  Please refer to
           the specific routine for a list of valid options.

        :Raises:
         - *ValueError* If the *file_format* requested does not match any of
           the available supported formats a *ValueError* is raised.
        """

        if file_format.lower() not in self._supported_formats.keys():
            raise ValueError("File format %s not available." % file_format)

        getattr(self, 'write_%s' % file_format.lower())(path, **kwargs)

    def write_geoclaw(self, path, force=False, skip=True, verbose=False,
                      fill_dict={}, **kwargs):
        r"""Write out a GeoClaw formatted storm file

        GeoClaw storm files are read in by the GeoClaw Fortran code.

        :Input:
         - *path* (string) Path to the file to be written.
         - *skip* (bool) Skip a time if NaNs are found and are not replaced.  
            Default is `True`.
         - *force* (bool) Force output of storm even if there is missing data.
            Default is `False`.
         - *verbose* (bool) Print out additional information when writing.
            Default is `False`.
         - *fill_dict* (dict) Dictionary of functions to use to fill in missing
            data represented by NaNs.  The keys are the field to be filled and
            the function signature should be `my_func(t, storm)` where t is the
            time of the forecast and `storm` is the storm object.  If the
            field remains a NaN or a function is not provided these lines will
            be assumed redundant and will be ommitted.  Note that the older 
            keyword arguments are put in this dictionary.  Currently the one
            default function is for `storm_radius`, which sets the value to 
            500 km.
        """

        # If a filling function is not provided we will provide some defaults
        fill_dict.update({"storm_radius": lambda t, storm: 500e3})
        # Handle older interface that had specific fill functions
        if "max_wind_radius_fill" in kwargs.keys():
            fill_dict.update(
                {"max_wind_radius": kwargs['max_wind_radius_fill']})
        if "storm_radius_fill" in kwargs.keys():
            fill_dict.update({"storm_radius": kwargs['storm_radius_fill']})

        # Loop through each line of data and if the line is valid, perform the
        # necessary work to write it out.  Otherwise either raise an exception
        # or skip it
        num_casts = 0
        data = []
        for n in range(len(self.t)):
            if self.t[n] == self.t[n - 1]:
                # Skip this time
                continue

            # Check each value we need for this time to make sure it is valid
            valid = True
            for name in ["max_wind_speed", "central_pressure",
                         "max_wind_radius", "storm_radius"]:
                if np.isnan(getattr(self, name)[n]):
                    if name in fill_dict.keys():
                        # Fill value with function provided
                        getattr(self, name)[n] = fill_dict[name](
                            self.t[n], self)
                    elif skip:
                        # Skip this line
                        valid = False
                        if verbose:
                            # Just warn that a NaN was found but continue
                            msg = ("*** WARNING:  The value {} at {} is a " +
                                   "NaN. Skipping this line.")
                            warnings.warn(msg.format(name, self.t[n]))
                    elif not force:
                        # If we are not asked to force to write raise an
                        # exception given the NaN
                        msg = ("The value {} at {} is a NaN and the storm " +
                               "will not be written in GeoClaw format.  If " +
                               "you want to fill in the value provide a " +
                               "function or set `force=True`.")
                        raise ValueError(msg.format(name, self.t[n]))
            if not valid:
                continue

            # Succeeded, add this time to the output
            num_casts += 1
            data.append(np.empty(7))

            # If we do not have a time offset use the first valid row as the
            # offset time
            if self.time_offset is None:
                self.time_offset = self.t[n]

            # Time
            if not isinstance(self.time_offset, float):
                data[-1][0] = (self.t[n] - self.time_offset).total_seconds()
            else:
                data[-1][0] = self.t[n] - self.time_offset
            # Eye-location
            data[-1][1:3] = self.eye_location[n, :]
            # Max wind speed
            data[-1][3] = self.max_wind_speed[n]
            # Max wind radius
            data[-1][4] = self.max_wind_radius[n]
            # Central pressure
            data[-1][5] = self.central_pressure[n]
            # Outer storm radius
            data[-1][6] = self.storm_radius[n]

        # Write out file
        format_string = ("{:19,.8e} " * 7)[:-1] + "\n"
        try:
            with open(path, "w") as data_file:
                # Write header
                data_file.write(f"{num_casts}\n")
                if isinstance(self.time_offset, np.datetime64):
                    data_file.write(\
                      f"{np.datetime_as_string(self.time_offset,unit='s')}\n\n")
                else:
                    data_file.write(f"{str(self.time_offset)}\n\n")

                # Write data lines
                for line in data:
                    data_file.write(format_string.format(*line))

        except Exception as e:
            # If an exception occurs clean up a partially generated file
            Path.unlink(path, missing_ok=True)
            raise e

    def write_atcf(self, path, verbose=False):
        r"""Write out a ATCF formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out ATCF files is not implemented ",
                                   "yet but is planned for a future release."))
        try:
            with open(path, 'w') as data_file:
                for n in range(len(self.t)):
                    data_file.write("".join((", " * 2,
                                             "%s" % seconds2date(self.t[n]),
                                             ", " * 4,
                                             "%s" % (int(self.eye_location[n, 0] *
                                                         10.0)),
                                             ", ",
                                             "%s" % (int(self.eye_location[n, 1] *
                                                         10.0)),
                                             ", ",
                                             "%s" % self.max_wind_speed[n],
                                             ", ",
                                             "%s" % self.central_pressure[n],
                                             ", ",
                                             ", " * 8,
                                             "%s" % self.storm_radius[n],
                                             ", ",
                                             "%s" % self.max_wind_radius[n],
                                             ", " * 10,
                                             "\n")))
        except Exception as e:
            # Remove possiblly partially generated file if not successful
            Path.unlink(path, missing_ok=True)
            raise e

    def write_hurdat(self, path, verbose=False):
        r"""Write out a HURDAT formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out hurdat files is not ",
                                   "implemented yet but is planned for a ",
                                   "future release."))
        try:
            with open(path, 'w') as data_file:
                data_file.write('%s %s %s' % ("Date", "Hurricane Name",
                                              "Indicator"))
                for n in range(self.t.shape[0]):

                    latitude = float(self.eye_location[n, 0])
                    longitude = float(self.eye_location[n, 1])

                    # Convert latitude to proper Hurdat format e.g 12.0N
                    if latitude > 0:
                        latitude = str(np.abs(latitude)) + 'N'
                    else:
                        latitude = str(np.abs(latitude)) + 'S'

                    # Convert longitude to proper Hurdat format e.g 12.0W
                    if longitude > 0:
                        longitude = str(np.abs(longitude)) + 'E'
                    else:
                        longitude = str(np.abs(longitude)) + 'W'

                    data_file.write("".join(("%s" % self.seconds2date(
                        self.t[n])[0:-2],
                        "%s00" % self.seconds2date(
                        self.t[n])[-2:],
                        ", " * 3,
                        "%s" % (latitude),
                        ", ",
                        "%s" % (longitude),
                        ", ",
                        "%s" % self.max_wind_speed[n],
                        ", ",
                        "%s" % self.central_pressure[n],
                        ", ",
                        "%s" % self.storm_radius[n],
                        ", ",
                        "%s" % self.max_wind_radius[n],
                        ", " * 10,
                        "\n")))
        except Exception as e:
            # Remove possiblly partially generated file if not successful
            Path.unlink(path, missing_ok=True)
            raise e

    def write_jma(self, path, verbose=False):
        r"""Write out a JMA formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out JMA files is not implemented ",
                                   "yet but is planned for a future release."))
        try:
            with open(path, 'w') as data_file:
                for n in range(self.t.shape[0]):
                    data_file.write("".join(("%s" % self.seconds2date(self.t[n]),
                                             " " * 4,
                                             "%s" % (int(self.eye_location[n, 0] *
                                                         10.0)),
                                             ", ",
                                             "%s" % (int(self.eye_location[n, 1] *
                                                         10.0)),
                                             ", ",
                                             "%s" % self.max_wind_speed[n],
                                             ", ",
                                             "%s" % self.central_pressure[n],
                                             ", ",
                                             ", " * 8,
                                             "%s" % self.storm_radius[n],
                                             ", ",
                                             "%s" % self.max_wind_radius[n],
                                             ", " * 10,
                                             "\n")))
        except Exception as e:
            # Remove possiblly partially generated file if not successful
            Path.unlink(path, missing_ok=True)
            raise e

    def write_imd(self, path, verbose=False):
        r"""Write out an IMD formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
        """
        raise NotImplementedError(("Writing out IMD files is not implemented ",
                                   "yet but is planned for a future release."))

    def write_tcvitals(self, path, verbose=False):
        r"""Write out an TCVITALS formatted storm file

        :Input:
         - *path* (string) Path to the file to be written.
         - *verbose* (bool) Print out additional information when writing.
         """

        raise NotImplementedError(("Writing in TCVITALS files is not",
                                   "implemented yet but is planned for a ",
                                   "future release."))

    def write_data(self, path, dim_mapping=None, var_mapping=None, verbose=False):
        r"""
         """
        
        # Only one format right now
        _data_file_format_mapping = {'ascii': 1, 'nws12': 1, "owi": 1,
                                     'netcdf': 2, 'nws13': 2}

        _window_type_mapping = {None: 0, 'none': 0, 'custom': 1}

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

        if isinstance(self.window_type, int):
            window_type = self.window_type
        elif isinstance(self.window_type, str):
            if (self.window_type.lower() in 
                                    _window_type_mapping.keys()):
                window_type = _window_type_mapping[
                                          self.window_type.lower()]
            else:
                raise TypeError(f"Unknown window type" +
                                f" '{self.window_type}' provided.")
        else:
            raise TypeError(f"Unknown window type" +
                            f" '{self.window_type}' provided.")

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
            data_file.write(f"{str(window_type).ljust(20)} # Window type\n")
            data_file.write(f"{str(self.ramp_width).ljust(20)} # Ramp width\n")
            if window_type > 0:
                data_file.write(f"{str(self.window)[1:-1].ljust(20)} # Window\n")
            else:
                data_file.write(f"{str(None).ljust(20)} # Window\n")
            data_file.write("\n")
            data_file.write("# Format Data Information\n")
            if file_format == 1:
                # Check number of file paths
                if len(self.file_paths)%2 != 0:
                        raise ValueError("The number of files should be even, " + 
                                         "one for pressure and wind, for each " +
                                         "resolution provided.")
                data_file.write("\n")
            elif file_format == 2:
                # Get dimension mapping
                _dim_mapping = util.get_netcdf_names(self.file_paths[0], 
                                                     lookup_type='dim',
                                                     user_mapping=dim_mapping,
                                                     verbose=verbose)

                # Get variable mapping
                _var_mapping = util.get_netcdf_names(self.file_paths[0], 
                                                     lookup_type='var',
                                                     user_mapping=var_mapping,
                                                     verbose=verbose)
                data_file.write(f"{str(_dim_mapping['x'])} ")
                data_file.write(f"{str(_dim_mapping['y'])} ")
                data_file.write(f"{str(_dim_mapping['t'])}\n")
                data_file.write(f"{str(_var_mapping['wind_u'])} ")
                data_file.write(f"{str(_var_mapping['wind_v'])} ")
                data_file.write(f"{str(_var_mapping['pressure'])}\n")
                data_file.write("\n")

                if len(self.file_paths) != 1:
                    raise ValueError(f"Expected 1 path for NetCDF format, " + 
                                     f"got {len(self.file_paths)}")

            # Write paths
            data_file.write("# File paths\n")
            for path in self.file_paths:
                data_file.write(f"{path}\n")


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
                        (speeds < 30) * -1 +
                        (speeds >= 64) * (speeds < 83) * 1 +
                        (speeds >= 83) * (speeds < 96) * 2 +
                        (speeds >= 96) * (speeds < 113) * 3 +
                        (speeds >= 113) * (speeds < 135) * 4 +
                        (speeds >= 135) * 5)
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
# Model field construction - Models supported are
#  - Holland 1980 ('HOLLAND_1980') [1]
#  - Holland 2010 ('HOLLAND_2010') [2]
#  - Chavas, Lin, Emmanuel ('CLE_2015') [3]
# *TODO* - Add citations

# Dictionary of models.  Keys are function names, values are the proper name
# and a citation to the model
_supported_models = {"holland_1980": ["Holland 1980", "Holland, G. J. An Analytic Model of the Wind and Pressure Profiles in Hurricanes. Monthly Weather Review 108, 1212-1218 (1980)."],
                     "holland_2010": ["Holland 2010", "Holland, G. J., Belanger, J. I. & Fritz, A. A Revised Model for Radial Profiles of Hurricane Winds. Monthly Weather Review 138, 4393-4393 (2010)."],
                     "cle_2015": ["Chavas, Lin, Emmanuel 2015", "Chavas, D. R., Lin, N. & Emanuel, K. A Model for the Complete Radial Structure of the Tropical Cyclone Wind Field. Part I: Comparison with Observed Structure*. https://doi.org.ezproxy.cul.columbia.edu/10.1175/JAS-D-15-0014.1 72, 3647-3662 (2015)."]}


# In the case where the field is not rotationally symmetric then the r value
# defines the x and y axis extents.
def construct_fields(storm, r, t, model="holland_1980"):
    r""""""

    if model.lower() not in _supported_models.keys():
        raise ValueError("Model %s not available." % model)

    return getattr(sys.modules[__name__], model.lower())(storm, x, t)


# Specific implementations
def holland_1980(storm, r, t):
    r""""""
    raise NotImplementedError("Holland 1980 model has not been implemeted.")
    return None, None


def holland_2010(storm, r, t):
    r""""""
    raise NotImplementedError("Holland 2010 model has not been implemeted.")
    return None, None


def cle_2015(storm, r, t):
    r""""""
    raise NotImplementedError("CLE 2015 model has not been implemeted.")
    return None, None


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
                           coords={'t': getattr(storm_fill, 't')},
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
                           coords={'t': getattr(storm_targ, 't')},
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


# =============================================================================
# Utility functions
def available_formats():
    r"""Construct a string suitable for listing available storm file formats.
    """
    output = "Available Formats: (Function, Name, Citation)\n"
    for (model, values) in Storm._supported_formats.items():
        output = "".join((output, "%s: %s %s\n" % (values[0], model,
                                                   values[1])))
    return output


def available_models():
    r"""Construct a string suitable for listing available storm models.
    """
    output = "Function, Name, Citation\n"
    for (model, values) in _supported_models.items():
        output = "".join((output, "%s: %s %s\n" % (values[0], model,
                                                   values[1])))
    return output


def make_multi_structure(path):
    r"""Create a dictionary of Storm objects for ATCF files with multiple storm tracks in them
    """
    with open(path, 'r') as f:
        lines = f.readlines()
        curTime = "test"
        curTrack = "test"
        os.mkdir("Clipped_ATCFs")
        stormDict = {}
        for line in lines:
            lineArr = line.split(", ")
            if curTime in lineArr[2]:
                if curTrack in lineArr[4]:
                    fileWrite.writelines(line)
                else:
                    fileWrite.close()
                    stormDict[curTime].update({curTrack: Storm(path=os.path.join(os.path.expandvars(
                        os.getcwd()), "Clipped_ATCFs", curTime, curTrack), file_format="ATCF")})
                    curTrack = lineArr[4]
                    fileWrite = open("Clipped_ATCFs/" +
                                     curTime + "/" + curTrack, 'w')
                    fileWrite.writelines(line)
            else:
                if curTime != "test":
                    fileWrite.close()
                    stormDict[curTime].update({curTrack: Storm(path=os.path.join(os.path.expandvars(
                        os.getcwd()), "Clipped_ATCFs", curTime, curTrack), file_format="ATCF")})
                curTime = lineArr[2]
                curTrack = lineArr[4]
                stormDict[curTime] = {}
                os.mkdir("Clipped_ATCFs/" + curTime)
                fileWrite = open("Clipped_ATCFs/" +
                                 curTime + "/" + curTrack, 'w')
                fileWrite.writelines(line)
    return stormDict


class DataDerivedStorms(object):
    """
    Storm object to accommodate storms in the Oceanweather Inc. format.
    The data is structured as wind and pressure files that contain the wind
    and pressure fields for time steps every XX minutes.

    This storm object is initialized to read in both files and parse the data
    into a netcdf file for reading with the data_storm_module.f90

      :Attributes:
     - *t* (list(np.datetime64)) Contains the time at which each entry of
       the other arrays are at.  These are expected to be *datetime* objects.
       Note that when written some formats require a *time_offset* to be set.
     - *wind_speed* (ndarray(:, :)) Wind speeds defined in every record, such
       as 34kt, 50kt, 64kt, etc and their radii. Default units are meters/second
       and meters.
     - *pressure* (ndarray(:,,:)) Pressure arrays every 15 minutes for the region
        of interest, Default units are bars
     - *u* (ndarray(:,:)) Wind velocity in the x direction in arrays for every 15 minutes
        Default units are m/s
     - *v* (ndarray(:,:)) Wind velocity in the y direction in arrays for every 15 minutes
        Default units are m/s
    - *latitude* (ndarray(:,:)) Array of latitudes for the wind and pressure arrays
    - *longitude* (ndarray(:,:)) Array of longitudes for the wind and pressure arrays

    :Initialization:
     1. Read in existing wind and pressure files at *path*.
     2. Construct a data derived storm object and parse data into attributes using
        data_storms.py for the parsing functions

    :Input:
     - *path* (string) Path to file to be read in if requested.
     - *kwargs* (dict) Other key-word arguments are passed to the appropriate
       read routine.
    """

    def __init__(self, path, wind_file_ext='WND', pressure_file_ext='PRE'):
        """
        This routine creates the DataStorm object, and loads the data from
        file using data_storms.py

        :param path: location and name of wind and pressure files
        :param wind_file_ext: Extension for the wind file, either 'WND' or 'WIN'
        :param pressure_file_ext: Extension for the pressure file "PRE" or 'pre'
        """
        # Set wind and pressure file extensions
        self.wind_ext = wind_file_ext
        self.pres_ext = pressure_file_ext

        # Read in wind and pressure data from original file formats
        self.wind_data = data_storms.read_oceanweather(path, self.wind_ext)
        self.pressure_data = data_storms.read_oceanweather(path, self.pres_ext)

        # Initialize instance variables
        self.t = None # Placeholder for array of time steps included in original data
        self.lat = None # Placeholder for array of latitudes of the storm domain
        self.lon = None # Placeholder for array of longitudes of the storm domain
        self.u = [] # Placeholder for wind data in x direction
        self.v = [] # Placeholder for wind data in y direction
        self.pressure = [] # Placeholder for pressure data
        self.wind_speed = [] # Placeholder for wind speed

    def parse_data(self,landfall_time):
        
        import numpy as np
        """
        This method processes the wind and pressure data matrices, extracts the data for each time step,
        and stores the processed data in instance variables.

        It uses functions from the data_storms module to extract time steps, latitude, and longitude arrays,
        as well as to process wind and pressure matrices into 2D arrays.

        The processed data is then appended to instance variables u, v, pressure, and wind_speed.
    
        :return: None
        """
        # Extract time and coordinate arrays
        self.t = data_storms.time_steps(self.wind_data, landfall_time)
        self.lat, self.lon = data_storms.get_coordinate_arrays(self.wind_data[0])

        # Iterate over all data, each index = an individual time step
        for winds, pressures in zip(self.wind_data, self.pressure_data):
            # Process the wind and pressure matrix per time index into 2d arrays
            u = data_storms.process_data(winds, start_idx=0)
            v = data_storms.process_data(winds, start_idx=(winds.iLat * winds.iLong))
            p = data_storms.process_data(pressures, start_idx=0)
            # Put into lists for easy writing to a xarray dataarray
            self.u.append(u)
            self.pressure.append(p)
            self.v.append(v)
            # Calculate wind speed from each directional component
            self.wind_speed.append(np.sqrt(u**2 + v**2))


    def write_data_derived(self, filename):
        """
        This method writes the derived wind, pressure, and wind speed data to a NetCDF file using xarray.

        :param filename: The name of the output NetCDF file (without the extension).
        :return: None
        """
        import xarray as xr

        windx = numpy.array(self.u)
        windy = numpy.array(self.v)
        pressure = numpy.array(self.pressure)*100 # Convert to mbars
        speed = numpy.array(self.wind_speed)
        time = numpy.array(self.t)

        # Create a dataset with xarray with derived data
        ds = xr.Dataset(data_vars={'u': (('time', 'lat', 'lon'), windx),
                                   'v': (('time', 'lat', 'lon'), windy),
                                   'speed': (('time', 'lat', 'lon'), speed),
                                   'pressure': (('time', 'lat', 'lon'), pressure),
                                   },
                        coords={'lat': self.lat,
                                'lon': self.lon,
                                'time': time})
        # Save the dataset to netcdf format
        ds.to_netcdf(filename + '.nc')











if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Positional argument
    parser.add_argument("path", help="Path to storm file to be read in")

    # Optional arguments
    parser.add_argument("-f", "--from", default="atcf", dest="input_format",
                        help="Format to convert from, defaults to 'atcf'")
    parser.add_argument("-o", "--output", default="geoclaw.storm",
                        dest="output_path",
                        help="Output path, default to 'geoclaw.storm'")
    parser.add_argument("-t", "--to", default="geoclaw",
                        dest="output_format",
                        help="Format to convert to, defaults to 'geoclaw'")
    parser.add_argument("-v", "--verbose",
                        help="Increase verbosity of output",
                        action="store_true")

    args = parser.parse_args()
    input_storm = Storm(args.path, file_format=args.input_format,
                        verbose=args.verbose)
    input_storm.write(args.output_path, file_format=args.output_format,
                      verbose=args.verbose)
