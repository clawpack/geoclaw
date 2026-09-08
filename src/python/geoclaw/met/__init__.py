# encoding: utf-8
r"""Meteorological (storm) forcing model for GeoClaw.

This package holds GeoClaw's storm/meteorological forcing object model: storm
tracks (:class:`~clawpack.geoclaw.met.track.StormTrack`), parametric wind/pressure
models (:class:`~clawpack.geoclaw.met.parametric.ParametricMetForcing`), gridded
file-backed forcing (:class:`~clawpack.geoclaw.met.gridded.GriddedMetForcing`),
and the historical :class:`~clawpack.geoclaw.met.storm.Storm` wrapper.

It was previously named :mod:`clawpack.geoclaw.surge`, which now remains as a
deprecated alias (importing it emits a :class:`DeprecationWarning`).  The
forcing-configuration data object lives in :mod:`clawpack.geoclaw.data` as
:class:`~clawpack.geoclaw.data.SurgeData`, aliased there (and re-exported here) as
``MetData``.  The on-disk ``surge.data`` filename and the Fortran modules are
unchanged.

The ``plot`` submodule (:mod:`clawpack.geoclaw.met.plot`) is imported lazily to
avoid a hard ``matplotlib`` dependency for non-plotting use.
"""

from clawpack.geoclaw.met.storm import (  # noqa: F401
    Storm, construct_fields, available_formats, available_models)
from clawpack.geoclaw.met.track import Track, StormTrack  # noqa: F401
from clawpack.geoclaw.met.track import (  # noqa: F401
    iter_hurdat, iter_ibtracs, catalog_hurdat, catalog_ibtracs)
from clawpack.geoclaw.met.tools import iter_atcf  # noqa: F401
from clawpack.geoclaw.met import reconstruction  # noqa: F401
from clawpack.geoclaw.met.parametric import ParametricMetForcing  # noqa: F401
from clawpack.geoclaw.met.gridded import GriddedMetForcing  # noqa: F401
from clawpack.geoclaw.met.data_storms import OWIData  # noqa: F401
from clawpack.geoclaw.data import SurgeData, MetData  # noqa: F401

__all__ = [
    'Storm', 'Track', 'StormTrack',
    'iter_atcf', 'iter_hurdat', 'iter_ibtracs',
    'catalog_hurdat', 'catalog_ibtracs',
    'reconstruction',
    'ParametricMetForcing', 'GriddedMetForcing', 'OWIData',
    'SurgeData', 'MetData',
    'construct_fields', 'available_formats', 'available_models',
    'plot',
]
