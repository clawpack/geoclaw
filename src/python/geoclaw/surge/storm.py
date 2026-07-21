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

:Object model:
    ``Storm`` is a thin backwards-compatibility wrapper over the meteorological
    forcing object model introduced by the refactor:
    :class:`~clawpack.geoclaw.surge.track.StormTrack`,
    :class:`~clawpack.geoclaw.surge.parametric.ParametricMetForcing`, and
    :class:`~clawpack.geoclaw.surge.gridded.GriddedMetForcing`.  Every legacy
    attribute is exposed as a property delegating to the underlying object.
"""

import sys
import argparse

import numpy as np

# Re-exported for backwards compatibility of the historical import paths and
# the module-level public surface.
from clawpack.geoclaw.surge.track import (  # noqa: F401
    Track, StormTrack,
    ATCF_basins, TCVitals_Basins, TC_designations, hurdat_special_entries,
    missing_data_warning_str, missing_necessary_data_warning_str,
    NoDataError, fill_rad_w_other_source)
from clawpack.geoclaw.surge.parametric import ParametricMetForcing
from clawpack.geoclaw.surge.gridded import GriddedMetForcing


def _track_delegate(name):
    r"""Property proxying attribute *name* to ``self.track``."""
    def getter(self):
        return getattr(self.track, name)

    def setter(self, value):
        setattr(self.track, name, value)

    return property(getter, setter)


def _gridded_delegate(name):
    r"""Property proxying attribute *name* to ``self.gridded``."""
    def getter(self):
        return getattr(self.gridded, name)

    def setter(self, value):
        setattr(self.gridded, name, value)

    return property(getter, setter)


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

    ``Storm`` is a thin compatibility wrapper: it holds a
    :class:`~clawpack.geoclaw.surge.track.StormTrack` (plus a
    :class:`~clawpack.geoclaw.surge.parametric.ParametricMetForcing` for the
    parameterized/track path) and a
    :class:`~clawpack.geoclaw.surge.gridded.GriddedMetForcing` for the gridded
    (data) path, and exposes every legacy attribute as a property delegating to
    the appropriate underlying object.

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

        # The wrapper always holds a track (parameterized/track path) and a
        # gridded forcing (data path).  Both share one bookkeeping object so
        # the shared fields (time_offset, file_paths, file_format) stay
        # consistent regardless of which path is active.
        self.track = StormTrack()
        self.parametric = ParametricMetForcing(track=self.track)
        self.gridded = GriddedMetForcing(meta=self.track.meta)

        if path is not None:
            self.read(path, file_format=file_format, **kwargs)

    # ==========================================================================
    #  Delegating attribute surface

    # Track fields
    t = _track_delegate("t")
    eye_location = _track_delegate("eye_location")
    max_wind_speed = _track_delegate("max_wind_speed")
    max_wind_radius = _track_delegate("max_wind_radius")
    central_pressure = _track_delegate("central_pressure")
    storm_radius = _track_delegate("storm_radius")
    wind_speeds = _track_delegate("wind_speeds")
    name = _track_delegate("name")
    basin = _track_delegate("basin")
    ID = _track_delegate("ID")
    classification = _track_delegate("classification")
    event = _track_delegate("event")

    # Shared bookkeeping (track and gridded share the same ``_Meta``).
    time_offset = _track_delegate("time_offset")
    file_paths = _track_delegate("file_paths")
    file_format = _track_delegate("file_format")

    # Gridded / forcing-control fields
    scaling = _gridded_delegate("scaling")
    storm_time_scale = _gridded_delegate("storm_time_scale")
    crop_extent = _gridded_delegate("crop_extent")
    ramp_width = _gridded_delegate("ramp_width")
    x_shift = _gridded_delegate("x_shift")
    y_shift = _gridded_delegate("y_shift")
    met_x_name = _gridded_delegate("met_x_name")
    met_y_name = _gridded_delegate("met_y_name")
    met_time_name = _gridded_delegate("met_time_name")
    met_lon_wrap = _gridded_delegate("met_lon_wrap")
    met_y_increasing = _gridded_delegate("met_y_increasing")
    met_fill_value = _gridded_delegate("met_fill_value")
    met_fill_action = _gridded_delegate("met_fill_action")
    met_time_offset = _gridded_delegate("met_time_offset")
    met_variable_map = _gridded_delegate("met_variable_map")

    # ==========================================================================
    #  Adoption helpers (keep the wrapper's object graph consistent post-read)
    def _adopt_track(self, track):
        r"""Adopt a ``StormTrack`` produced by a track-format reader."""
        self.track = track
        self.parametric.track = track
        self.gridded.meta = track.meta

    def _adopt_parametric(self, forcing):
        r"""Adopt a ``ParametricMetForcing`` produced by ``read_geoclaw``."""
        self.parametric = forcing
        self.track = forcing.track
        self.gridded.meta = forcing.track.meta

    def _adopt_gridded(self, forcing):
        r"""Adopt a ``GriddedMetForcing`` produced by ``read_data``."""
        self.gridded = forcing
        self.track.meta = forcing.meta
        self.parametric.track = self.track

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
        from pathlib import Path

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

    def read_geoclaw(self, path, **kwargs):
        r"""Read in a GeoClaw formatted storm file (delegates to
        :meth:`ParametricMetForcing.read_geoclaw`)."""
        self._adopt_parametric(ParametricMetForcing.read_geoclaw(path, **kwargs))

    def read_atcf(self, path, **kwargs):
        r"""Read in an ATCF formatted storm file (delegates to
        :meth:`StormTrack.read_atcf`)."""
        self._adopt_track(StormTrack.read_atcf(path, **kwargs))

    def read_hurdat(self, path, **kwargs):
        r"""Read in a HURDAT formatted storm file (delegates to
        :meth:`StormTrack.read_hurdat`)."""
        self._adopt_track(StormTrack.read_hurdat(path, **kwargs))

    def read_ibtracs(self, path, **kwargs):
        r"""Read in an IBTrACS formatted storm file (delegates to
        :meth:`StormTrack.read_ibtracs`)."""
        self._adopt_track(StormTrack.read_ibtracs(path, **kwargs))

    def read_jma(self, path, **kwargs):
        r"""Read in a JMA formatted storm file (delegates to
        :meth:`StormTrack.read_jma`)."""
        self._adopt_track(StormTrack.read_jma(path, **kwargs))

    def read_imd(self, path, **kwargs):
        r"""Read in an IMD formatted storm file (delegates to
        :meth:`StormTrack.read_imd`)."""
        self._adopt_track(StormTrack.read_imd(path, **kwargs))

    def read_tcvitals(self, path, **kwargs):
        r"""Read in a TCVITALS formatted storm file (delegates to
        :meth:`StormTrack.read_tcvitals`)."""
        self._adopt_track(StormTrack.read_tcvitals(path, **kwargs))

    def read_data(self, path, **kwargs):
        r"""Read in a data-derived (gridded) storm file (delegates to
        :meth:`GriddedMetForcing.read_data`)."""
        self._adopt_gridded(GriddedMetForcing.read_data(path, **kwargs))

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

    def write_geoclaw(self, path, **kwargs):
        r"""Write out a GeoClaw formatted storm file (delegates to
        :meth:`ParametricMetForcing.write_geoclaw`)."""
        self.parametric.write_geoclaw(path, **kwargs)

    def write_data(self, path, **kwargs):
        r"""Write out a data-derived (gridded) storm descriptor (delegates to
        :meth:`GriddedMetForcing.write_data`)."""
        self.gridded.write_data(path, **kwargs)

    # -- Track-format writers are not primary supported outputs (stubs) -------
    def write_atcf(self, path, verbose=False):
        r"""(Not implemented) Write out a ATCF formatted storm file."""
        raise NotImplementedError(("Writing out ATCF files is not implemented ",
                                   "yet but is planned for a future release."))

    def write_hurdat(self, path, verbose=False):
        r"""(Not implemented) Write out a HURDAT formatted storm file."""
        raise NotImplementedError(("Writing out hurdat files is not ",
                                   "implemented yet but is planned for a ",
                                   "future release."))

    def write_jma(self, path, verbose=False):
        r"""(Not implemented) Write out a JMA formatted storm file."""
        raise NotImplementedError(("Writing out JMA files is not implemented ",
                                   "yet but is planned for a future release."))

    def write_imd(self, path, verbose=False):
        r"""(Not implemented) Write out an IMD formatted storm file."""
        raise NotImplementedError(("Writing out IMD files is not implemented ",
                                   "yet but is planned for a future release."))

    def write_tcvitals(self, path, verbose=False):
        r"""(Not implemented) Write out a TCVITALS formatted storm file."""
        raise NotImplementedError(("Writing in TCVITALS files is not",
                                   "implemented yet but is planned for a ",
                                   "future release."))

    # =========================================================================
    #  Track plotting / categorization (delegated to the StormTrack)
    def plot(self, ax, *args, **kwargs):
        r"""Plot this storm's track (delegates to :meth:`StormTrack.plot`)."""
        return self.track.plot(ax, *args, **kwargs)

    def category(self, categorization="NHC", cat_names=False):
        r"""Categorize the storm (delegates to :meth:`StormTrack.category`)."""
        return self.track.category(categorization=categorization,
                                   cat_names=cat_names)


# =============================================================================
# Model field construction
#
# The Python parameterized-model routines below are stubs; the real
# implementations live in the GeoClaw Fortran code.  They are retained (and
# kept importable) but are no longer part of the advertised object surface.
# A future ``parametric`` evaluator will take a ``StormTrack`` and output
# wind/pressure on a space-time grid, mirroring the Fortran.

# Dictionary of Python model stubs.  Keys are function names, values are the
# proper name and a citation to the model.  Retained for ``construct_fields``.
_supported_models = {"holland_1980": ["Holland 1980", "Holland, G. J. An Analytic Model of the Wind and Pressure Profiles in Hurricanes. Monthly Weather Review 108, 1212-1218 (1980)."],
                     "holland_2010": ["Holland 2010", "Holland, G. J., Belanger, J. I. & Fritz, A. A Revised Model for Radial Profiles of Hurricane Winds. Monthly Weather Review 138, 4393-4393 (2010)."],
                     "cle_2015": ["Chavas, Lin, Emmanuel 2015", "Chavas, D. R., Lin, N. & Emanuel, K. A Model for the Complete Radial Structure of the Tropical Cyclone Wind Field. Part I: Comparison with Observed Structure*. https://doi.org.ezproxy.cul.columbia.edu/10.1175/JAS-D-15-0014.1 72, 3647-3662 (2015)."]}

# Parameterized-forcing models actually supported by the GeoClaw Fortran code
# (see met_forcing_refactor.md Section 7).  Labeled separately from any Python
# implementations (of which there are none yet).
_fortran_supported_models = {
    "holland80": "Holland 1980",
    "holland2008": "Holland 2008",
    "holland2010": "Holland 2010",
    "cle": "Chavas, Lin, Emanuel 2015",
    "slosh": "SLOSH",
    "rankine": "Rankine vortex",
    "modified_rankine": "Modified Rankine vortex",
    "demaria": "DeMaria",
    "willoughby": "Willoughby",
}

# Parameterized-forcing models implemented in Python.  None are functional yet;
# the Python Holland/CLE routines are stubs that raise ``NotImplementedError``.
_python_supported_models = {}


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
    r"""Construct a string listing the available parameterized-forcing models.

    The models are split into those supported by the GeoClaw Fortran code and
    those implemented in Python (currently none), both clearly labeled.
    """
    output = "Fortran-supported models: (Function, Name)\n"
    for (model, name) in _fortran_supported_models.items():
        output = "".join((output, "%s: %s\n" % (model, name)))
    output = "".join((output, "Python-supported models: (Function, Name)\n"))
    if _python_supported_models:
        for (model, name) in _python_supported_models.items():
            output = "".join((output, "%s: %s\n" % (model, name)))
    else:
        output = "".join((output, "(none implemented yet)\n"))
    return output


# ``make_multi_structure`` moved to the workflow-tools namespace; kept
# importable here for backwards compatibility.
from clawpack.geoclaw.surge.tools import make_multi_structure  # noqa: E402,F401


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
