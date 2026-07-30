# encoding: utf-8
r"""Deprecated alias for :mod:`clawpack.geoclaw.met`.

``clawpack.geoclaw.surge`` was renamed to :mod:`clawpack.geoclaw.met`.  Importing
anything from ``surge`` still works but emits a :class:`DeprecationWarning`; update
imports to ``clawpack.geoclaw.met`` and use
:class:`clawpack.geoclaw.data.MetData` in place of ``SurgeData``.  The on-disk
``surge.data`` filename and the Fortran modules are unchanged.

The warning is emitted here, at the package level, so it fires once for any
``clawpack.geoclaw.surge.*`` import (importing a submodule imports this package
first); the submodule shims themselves are silent.
"""
import warnings as _warnings

_warnings.warn(
    "clawpack.geoclaw.surge has been renamed to clawpack.geoclaw.met; the "
    "'surge' name is deprecated and will be removed in a future release. "
    "Update imports to 'clawpack.geoclaw.met'.",
    DeprecationWarning, stacklevel=2)

from clawpack.geoclaw import met as _met  # noqa: E402
from clawpack.geoclaw.met import (  # noqa: E402,F401
    Storm, Track, StormTrack, ParametricMetForcing, GriddedMetForcing,
    OWIData, SurgeData, MetData,
    construct_fields, available_formats, available_models)

__all__ = [
    'Storm', 'Track', 'StormTrack',
    'ParametricMetForcing', 'GriddedMetForcing', 'OWIData',
    'SurgeData', 'MetData',
    'construct_fields', 'available_formats', 'available_models',
    'plot',
]


def __getattr__(name):
    # Delegate any other public name to the canonical met package.
    return getattr(_met, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_met)))
