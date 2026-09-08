# encoding: utf-8
r"""Deprecated alias for :mod:`clawpack.geoclaw.met.reconstruction`.

Import from :mod:`clawpack.geoclaw.met.reconstruction` instead.  ``clawpack.geoclaw.surge``
is deprecated; a :class:`DeprecationWarning` is emitted once from
:mod:`clawpack.geoclaw.surge` when the package is imported.
"""
from clawpack.geoclaw.met import reconstruction as _target


def __getattr__(name):
    return getattr(_target, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_target)))
