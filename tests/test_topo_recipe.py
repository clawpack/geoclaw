#!/usr/bin/env python
# encoding: utf-8
"""Per-file input recipes: ``TopoFile`` and ``DTopoFile`` (clawpack/geoclaw#768).

A recipe is a file plus the settings ``setrun.py`` uses to turn it into what
Fortran loads.  ``Topography`` and ``DTopography`` hold data and I/O only.
Design: dev/design/topo_recipe.md.
"""

import warnings

import numpy as np
import pytest

from clawpack.geoclaw import dtopotools, topotools
from clawpack.geoclaw.data import DTopoData, TopographyData

pytestmark = pytest.mark.python

RECIPE = ("crop_extent", "coarsen", "buffer", "align",
          "x_shift", "y_shift", "z_shift", "negate_z")
DTOPO_RECIPE = ("x_shift", "y_shift", "z_shift", "negate_z")

FILE_X = (-100.0, -80.0)
FILE_Y = (20.0, 30.0)
DX = 0.125
CROP = [-95.0, -90.0, 22.0, 25.0]


def _field(X, Y):
    return -100.0 + 2.0 * (X - FILE_X[0]) + 3.0 * (Y - FILE_Y[0])


def _write_topo(tmp_path, topo_type):
    x = np.arange(FILE_X[0], FILE_X[1] + DX / 2, DX)
    y = np.arange(FILE_Y[0], FILE_Y[1] + DX / 2, DX)
    X, Y = np.meshgrid(x, y)
    Z = _field(X, Y)
    if topo_type == 4:
        xr = pytest.importorskip("xarray")
        pytest.importorskip("netCDF4")
        path = tmp_path / "topo.nc"
        xr.Dataset(
            {"z": (("lat", "lon"), Z, {"units": "m"})},
            coords={
                "lon": ("lon", x, {"units": "degrees_east",
                                   "standard_name": "longitude"}),
                "lat": ("lat", y, {"units": "degrees_north",
                                   "standard_name": "latitude"}),
            }).to_netcdf(path)
        return path
    t = topotools.Topography()
    t.x, t.y, t.Z = x, y, Z
    path = tmp_path / f"topo.tt{topo_type}"
    t.write(str(path), topo_type=topo_type)
    return path


def _write_dtopo(tmp_path):
    x = np.linspace(-99.0, -81.0, 37)
    y = np.linspace(21.0, 29.0, 17)
    X, Y = np.meshgrid(x, y)
    d = dtopotools.DTopography()
    d.x, d.y, d.X, d.Y = x, y, X, Y
    d.times = np.array([0.0, 1.0])
    d.dZ = np.array([0.0 * X, 0.1 * np.exp(-((X + 90) ** 2 + (Y - 25) ** 2))])
    path = tmp_path / "dtopo.dtt3"
    d.write(str(path), dtopo_type=3)
    return path


# --------------------------------------------------------------------------
# TopoFile and DTopoFile: same fields, defaults, typos, validation
# --------------------------------------------------------------------------

def _topofile():
    return topotools.TopoFile, "topo_type"


def _dtopofile():
    return dtopotools.DTopoFile, "dtopo_type"


RECIPE_CLASSES = [pytest.param(_topofile, id="TopoFile"),
                  pytest.param(_dtopofile, id="DTopoFile")]

KEYWORDS = dict(crop_extent=CROP, coarsen=2, buffer=1, align=[-95.0, 22.0],
                x_shift=1.0, y_shift=-1.0, z_shift=0.5, negate_z=True)


@pytest.mark.parametrize("which", RECIPE_CLASSES)
def test_recipe_defaults(which):
    cls, type_name = which()
    e = cls()
    assert e.path is None
    assert getattr(e, type_name) is None
    assert e.crop_extent is None
    assert e.coarsen == 1
    assert e.buffer == 0
    assert e.align is None
    assert (e.x_shift, e.y_shift, e.z_shift) == (0.0, 0.0, 0.0)
    assert e.negate_z is False
    assert e.datum is None


@pytest.mark.parametrize("which", RECIPE_CLASSES)
def test_recipe_positional_and_keywords(which):
    cls, type_name = which()
    e = cls("a.nc", 4, **KEYWORDS)
    assert e.path == "a.nc"
    assert getattr(e, type_name) == 4
    assert list(e.crop_extent) == CROP
    assert (e.coarsen, e.buffer) == (2, 1)
    assert list(e.align) == [-95.0, 22.0]
    assert (e.x_shift, e.y_shift, e.z_shift, e.negate_z) == (1.0, -1.0, 0.5, True)


@pytest.mark.parametrize("which", RECIPE_CLASSES)
def test_recipe_keywords_equal_attributes_after_creation(which):
    """Randy (2026-11-24): an entry can be built with keywords or by setting
    attributes after creation, as with ``Topography``; both give the same
    entry."""
    cls, type_name = which()
    by_keyword = cls(path="a.nc", **{type_name: 4}, **KEYWORDS)
    by_attribute = cls()
    by_attribute.path = "a.nc"
    setattr(by_attribute, type_name, 4)
    for name, value in KEYWORDS.items():
        setattr(by_attribute, name, value)
    assert by_keyword == by_attribute


@pytest.mark.parametrize("which", RECIPE_CLASSES)
def test_recipe_misspelled_attribute_raises(which):
    """Randy's #765 repro: ``crop_exent`` was silently a new attribute, so the
    file was written with no crop and Fortran loaded all of it."""
    cls, _ = which()
    e = cls("a.tt3", 3)
    with pytest.raises(AttributeError):
        e.crop_exent = CROP


@pytest.mark.parametrize("which", RECIPE_CLASSES)
def test_recipe_misspelled_keyword_raises(which):
    cls, _ = which()
    with pytest.raises(TypeError):
        cls("a.tt3", 3, crop_exent=CROP)


@pytest.mark.parametrize("which", RECIPE_CLASSES)
@pytest.mark.parametrize("name, value", [
    ("crop_extent", [1.0, 2.0, 3.0]),          # wrong length
    ("crop_extent", [2.0, 1.0, 3.0, 4.0]),     # x1 >= x2
    ("crop_extent", [1.0, 2.0, 4.0, 3.0]),     # y1 >= y2
    ("coarsen", 0),
    ("buffer", -1),
    ("buffer", 1.5),                           # a count of grid points
    ("align", [1.0, 2.0, 3.0]),
])
def test_recipe_rejects_invalid_values(which, name, value):
    cls, _ = which()
    with pytest.raises(ValueError):
        cls("a.tt3", 3, **{name: value})
    e = cls("a.tt3", 3)
    with pytest.raises(ValueError):
        setattr(e, name, value)


# --------------------------------------------------------------------------
# Topography / DTopography are data + I/O only
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name", RECIPE)
def test_topography_has_no_recipe_attribute(name):
    assert not hasattr(topotools.Topography(), name)


@pytest.mark.parametrize("name", RECIPE)
def test_dtopography_has_no_recipe_attribute(name):
    assert not hasattr(dtopotools.DTopography(), name)


@pytest.mark.parametrize("cls, recipe", [
    (topotools.Topography, "TopoFile"),
    (dtopotools.DTopography, "DTopoFile"),
])
@pytest.mark.parametrize("name", RECIPE)
def test_data_classes_refuse_setting_names(cls, recipe, name):
    """Code written against master before #768 set the settings on the data
    object.  As plain attributes they would be stored and silently ignored by
    read(), so they are refused, and the message says where they went."""
    with pytest.raises(AttributeError, match=recipe):
        setattr(cls(), name, 1)


# --------------------------------------------------------------------------
# TopoFile.read() / read_header()
# --------------------------------------------------------------------------

@pytest.mark.parametrize("topo_type", [2, 3, 4])
def test_topofile_read_matches_topography_read_keywords(tmp_path, topo_type):
    """``TopoFile.read()`` is ``Topography.read`` with the recipe passed as
    one-off keywords; the result carries no recipe."""
    path = _write_topo(tmp_path, topo_type)
    recipe = dict(crop_extent=CROP, coarsen=2, buffer=1,
                  align=[-94.875, 22.125], x_shift=2.0, y_shift=1.0,
                  z_shift=0.5, negate_z=True)
    got = topotools.TopoFile(str(path), topo_type, **recipe).read()
    want = topotools.Topography()
    want.read(str(path), topo_type=topo_type, **recipe)
    assert isinstance(got, topotools.Topography)
    np.testing.assert_array_equal(got.x, want.x)
    np.testing.assert_array_equal(got.y, want.y)
    np.testing.assert_array_equal(got.Z, want.Z)


@pytest.mark.parametrize("topo_type", [2, 3, 4])
def test_topography_read_shift_keywords(tmp_path, topo_type):
    """Shifts and ``negate_z`` act on one read and are not stored."""
    path = _write_topo(tmp_path, topo_type)
    plain = topotools.Topography()
    plain.read(str(path), topo_type=topo_type)
    t = topotools.Topography()
    t.read(str(path), topo_type=topo_type, x_shift=2.0, y_shift=1.0,
           z_shift=0.5, negate_z=True)
    np.testing.assert_array_equal(t.x, plain.x + 2.0)
    np.testing.assert_array_equal(t.y, plain.y + 1.0)
    np.testing.assert_allclose(t.Z, -plain.Z + 0.5, rtol=0, atol=1e-12)


@pytest.mark.parametrize("topo_type", [2, 3, 4])
def test_topofile_read_header_is_the_file(tmp_path, topo_type):
    """``read_header`` reports the file as stored: no crop, no shift."""
    path = _write_topo(tmp_path, topo_type)
    e = topotools.TopoFile(str(path), topo_type, crop_extent=CROP,
                           x_shift=2.0, y_shift=1.0)
    h = e.read_header()
    np.testing.assert_allclose(h.extent, [*FILE_X, *FILE_Y], rtol=0,
                               atol=1e-9)


# --------------------------------------------------------------------------
# DTopoFile.read() / read_header(), and the write guard until #15
# --------------------------------------------------------------------------

DTOPO_RECIPE_KW = dict(crop_extent=[-95.0, -85.0, 22.0, 27.0], coarsen=2,
                       buffer=1, align=[-94.5, 22.5], x_shift=2.0,
                       y_shift=1.0, z_shift=0.5, negate_z=True)


def test_dtopofile_read_matches_dtopography_read_keywords(tmp_path):
    path = _write_dtopo(tmp_path)
    got = dtopotools.DTopoFile(str(path), 3, **DTOPO_RECIPE_KW).read()
    want = dtopotools.DTopography()
    want.read(str(path), dtopo_type=3, **DTOPO_RECIPE_KW)
    assert isinstance(got, dtopotools.DTopography)
    np.testing.assert_array_equal(got.x, want.x)
    np.testing.assert_array_equal(got.y, want.y)
    np.testing.assert_array_equal(got.dZ, want.dZ)


def test_dtopography_read_keywords(tmp_path):
    """The settings act on one read and are not stored."""
    path = _write_dtopo(tmp_path)
    plain = dtopotools.DTopography()
    plain.read(str(path), dtopo_type=3)
    d = dtopotools.DTopography()
    d.read(str(path), dtopo_type=3, x_shift=2.0, y_shift=1.0)
    np.testing.assert_array_equal(d.x, plain.x + 2.0)
    np.testing.assert_array_equal(d.y, plain.y + 1.0)
    c = dtopotools.DTopography()
    c.read(str(path), dtopo_type=3, coarsen=2)
    np.testing.assert_array_equal(c.x, plain.x[::2])


def test_dtopofile_read_header_is_the_file(tmp_path):
    path = _write_dtopo(tmp_path)
    e = dtopotools.DTopoFile(str(path), 3, crop_extent=[-95.0, -85.0, 22.0, 27.0],
                             x_shift=2.0, y_shift=1.0)
    h = e.read_header()
    np.testing.assert_allclose(h.extent, [-99.0, -81.0, 21.0, 29.0],
                               rtol=0, atol=1e-9)


@pytest.mark.parametrize("name, value", [
    ("crop_extent", [-95.0, -85.0, 22.0, 27.0]),
    ("coarsen", 2),
    ("buffer", 1),
    ("align", [-94.5, 22.5]),
])
def test_dtopo_data_refuses_settings_fortran_lacks(tmp_path, name, value):
    """Fortran does not crop or coarsen dtopo yet (mandli/geoclaw#15), so
    ``make data`` stops instead of writing a file the run would reject."""
    path = _write_dtopo(tmp_path)
    with pytest.raises(NotImplementedError, match="#15"):
        _write_dtopo_data(tmp_path,
                          [dtopotools.DTopoFile(str(path), 3, **{name: value})])


# --------------------------------------------------------------------------
# What topofiles / dtopofiles accept
# --------------------------------------------------------------------------

def _write_topo_data(tmp_path, entries):
    td = TopographyData()
    td.topofiles = entries
    out = tmp_path / "topo.data"
    td.write(out_file=str(out), data_source="test")
    return out.read_text()


def test_topofiles_accept_topofile(tmp_path):
    path = _write_topo(tmp_path, 3)
    text = _write_topo_data(tmp_path, [topotools.TopoFile(str(path), 3,
                                                          crop_extent=CROP)])
    assert str(path) in text


def test_topofiles_list_form_is_not_deprecated(tmp_path):
    """``[topo_type, path]`` is the v5.14 form and every example uses it."""
    path = _write_topo(tmp_path, 3)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        list_text = _write_topo_data(tmp_path, [[3, str(path)]])
    assert list_text == _write_topo_data(
        tmp_path, [topotools.TopoFile(str(path), 3)])


@pytest.mark.parametrize("entry", ["topography", "dict"])
def test_topofiles_reject_unreleased_forms(tmp_path, entry):
    path = _write_topo(tmp_path, 3)
    if entry == "topography":
        e = topotools.Topography()
        e.path, e.topo_type = str(path), 3
    else:
        e = {"topo_path": str(path), "topo_type": 3}
    with pytest.raises(TypeError, match="TopoFile"):
        _write_topo_data(tmp_path, [e])


def _write_dtopo_data(tmp_path, entries):
    dd = DTopoData()
    dd.dtopofiles = entries
    out = tmp_path / "dtopo.data"
    dd.write(out_file=str(out), data_source="test")
    return out.read_text()


def test_dtopofiles_accept_dtopofile_and_list(tmp_path):
    path = _write_dtopo(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        list_text = _write_dtopo_data(tmp_path, [[3, str(path)]])
    assert list_text == _write_dtopo_data(
        tmp_path, [dtopotools.DTopoFile(str(path), 3)])


def test_dtopofiles_reject_dtopography(tmp_path):
    path = _write_dtopo(tmp_path)
    d = dtopotools.DTopography()
    d.path, d.dtopo_type = str(path), 3
    with pytest.raises(TypeError, match="DTopoFile"):
        _write_dtopo_data(tmp_path, [d])


def test_topofile_without_path_raises_at_write(tmp_path):
    with pytest.raises(ValueError, match="no path"):
        _write_topo_data(tmp_path, [topotools.TopoFile(topo_type=3)])


def test_dtopofile_without_path_raises_at_write(tmp_path):
    with pytest.raises(ValueError, match="no path"):
        _write_dtopo_data(tmp_path, [dtopotools.DTopoFile(dtopo_type=3)])
