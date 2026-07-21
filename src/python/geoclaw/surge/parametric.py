#!/usr/bin/env python

r"""Parameterized meteorological forcing.

``ParametricMetForcing`` holds a :class:`~clawpack.geoclaw.surge.track.Track`
(or :class:`~clawpack.geoclaw.surge.track.StormTrack`), a model identifier, and
model configuration.  It provides the GeoClaw parametric-forcing reader
(``read_geoclaw``) and writer (``write_geoclaw``), both moved verbatim from the
original ``Storm`` implementation so the emitted bytes stay identical.

``time_offset`` is Python-to-Fortran conversion metadata resolved at the write
boundary exactly as before (it is stored on the track's shared ``_Meta`` and
proxied here).
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from clawpack.geoclaw.surge.track import StormTrack


def _track_prop(name):
    r"""Build a property that proxies attribute *name* to ``self.track``."""
    def getter(self):
        return getattr(self.track, name)

    def setter(self, value):
        setattr(self.track, name, value)

    return property(getter, setter)


class ParametricMetForcing(object):
    r"""Forcing generated from a parameterized model referencing a track.

    :Attributes:
     - *track* (:class:`StormTrack`) The track/event the forcing references.
     - *model* (str) Model identifier (e.g. ``"holland80"``).  Not required by
       the GeoClaw writer, which emits the compact track-parameter file.
    """

    def __init__(self, track=None, model=None):
        self.track = track if track is not None else StormTrack()
        self.model = model

    # Track-field proxies so the verbatim reader/writer bodies (which reference
    # ``self.<field>``) operate on the underlying track.
    name = _track_prop("name")
    t = _track_prop("t")
    eye_location = _track_prop("eye_location")
    max_wind_speed = _track_prop("max_wind_speed")
    max_wind_radius = _track_prop("max_wind_radius")
    central_pressure = _track_prop("central_pressure")
    storm_radius = _track_prop("storm_radius")

    # Shared bookkeeping (proxied through the track's ``_Meta``).
    time_offset = _track_prop("time_offset")
    file_paths = _track_prop("file_paths")
    file_format = _track_prop("file_format")

    # =========================================================================
    # Read Routines
    @classmethod
    def read_geoclaw(cls, path, verbose=False):
        r"""Read in a GeoClaw formatted storm file

        GeoClaw storm files are read in by the Fortran code and are not meant
        to be human readable.

        :Input:
         - *path* (string) Path to the file to be read.
         - *verbose* (bool) Output more info regarding reading.
        """
        self = cls()

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
            # GeoClaw storm files write elapsed time as floating-point seconds,
            # so construct timedeltas with sub-second support rather than
            # passing floats directly to ``np.timedelta64``.
            time_deltas = pd.to_timedelta(data[:, 0], unit="s").to_numpy()
            self.t = self.time_offset + time_deltas
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

        return self

    # =========================================================================
    # Write Routines
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

        # Get around the mutable-default-argument problem for the fill_dict
        if fill_dict is None:
            fill_dict = {}
        else:
            fill_dict = dict(fill_dict)

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
            if n > 0 and self.t[n] == self.t[n - 1]:
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
                delta = self.t[n] - self.time_offset
                if hasattr(delta, "total_seconds"):
                    data[-1][0] = float(delta.total_seconds())
                else:
                    if hasattr(delta, "values"):
                        delta = delta.values
                    if isinstance(delta, np.ndarray):
                        delta = delta.item()
                    data[-1][0] = float(pd.to_timedelta(delta).total_seconds())
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
