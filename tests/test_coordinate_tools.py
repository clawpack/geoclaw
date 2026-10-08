#!/usr/bin/env python
# encoding: utf-8
"""Tests for the product-/format-neutral coordinate primitives."""

import itertools

import numpy
import pytest

from clawpack.geoclaw.coordinate_tools import (
    is_geographic_lon,
    classify_lon_axis,
    resolve_wrap,
    _compute_lon_entries,
    axis_file_slice,
    crop_indices,
    EmptyCropWindow,
)


def _ref_crop_window(coords, lo, hi, delta, coarsen, buffer, align):
    """Original ``Topography.crop`` index math, inlined verbatim as an oracle.

    Returns the half-open ``(lower, upper)`` for one axis, or ``None`` when the
    request does not overlap, matching the pre-refactor behavior that
    ``crop_indices`` must reproduce exactly.
    """
    delta_new = delta * coarsen
    try:
        lower = (coords >= lo).nonzero()[0][0]
        upper = (coords <= hi).nonzero()[0][-1]
    except IndexError:
        return None
    if (coarsen > 1) and (align is not None):
        vs = numpy.array([coords[lower + i] for i in range(coarsen)])
        offsets = (vs - align) / delta_new
        offsets_frac = offsets - numpy.round(offsets)
        ioffset = numpy.argmin(abs(offsets_frac))
        lower = lower + ioffset
        upper = upper - numpy.remainder(upper - lower, coarsen)
    lower = numpy.maximum(0, lower - buffer * coarsen)
    upper = numpy.minimum(len(coords) - 1, upper + buffer * coarsen) + 1
    return int(lower), int(upper)


@pytest.mark.python
def test_is_geographic_lon_degree_ranges():
    assert is_geographic_lon(-190.0, -130.0) is True   # extended-frame lon-lat
    assert is_geographic_lon(0.0, 360.0) is True
    assert is_geographic_lon(-180.0, 180.0) is True


@pytest.mark.python
def test_is_geographic_lon_projected():
    # Large (meters-like) coordinates -> projected, never wrap.
    assert is_geographic_lon(0.0, 500000.0) is False
    # Length units on the axis -> projected even for small coords.
    assert is_geographic_lon(0.0, 10.0, units="km") is False
    assert is_geographic_lon(0.0, 10.0, units="m") is False
    # Projection standard_name -> projected.
    assert is_geographic_lon(0.0, 10.0,
                             std_name="projection_x_coordinate") is False
    # Small unit-less coords default to geographic (heuristic can't tell).
    assert is_geographic_lon(0.0, 10.0) is True


@pytest.mark.python
def test_compute_lon_entries_single_and_split():
    # Non-crossing crop within a [-180,180] file -> one identity entry.
    entries = _compute_lon_entries(-180.0, 180.0, -170.0, -150.0)
    assert entries == [(-170.0, -150.0, 0.0)]

    # Crop crossing -180 -> two entries with offsets {0, -360}.
    entries = _compute_lon_entries(-180.0, 180.0, -190.0, -170.0, max_gap=1.0)
    offsets = sorted(o for *_, o in entries)
    assert offsets == [-360.0, 0.0]
    pairs = {(round(a, 6), round(b, 6), o) for a, b, o in entries}
    assert (-180.0, -170.0, 0.0) in pairs
    assert (170.0, 180.0, -360.0) in pairs


@pytest.mark.python
def test_compute_lon_entries_no_wrap_disables_split():
    # allow_wrap=False (projected axis): only the identity offset is tried, so a
    # crop reaching past the file cannot be covered -> ValueError.
    with pytest.raises(ValueError):
        _compute_lon_entries(-180.0, 180.0, -190.0, -170.0,
                             max_gap=1.0, allow_wrap=False)


@pytest.mark.python
def test_crop_indices_matches_original_crop_math():
    """crop_indices reproduces the pre-refactor Topography.crop index math.

    Behavior-preservation Vet: sweep a grid of crop bounds x coarsen x buffer x
    align against the inlined original algorithm; the windows must be identical.
    """
    delta = 0.25
    coords = numpy.arange(0.0, 20.0 + delta / 2, delta)   # 0.0 .. 20.0
    los = [0.0, 1.1, 3.0, 7.6, 12.5]
    his = [4.0, 9.9, 15.0, 20.0]
    coarsens = [1, 2, 3, 4]
    buffers = [0, 1, 3]
    aligns = [None, 0.0, 0.5, 1.0]

    for lo, hi, coarsen, buffer, align in itertools.product(
            los, his, coarsens, buffers, aligns):
        if hi <= lo:
            continue
        ref = _ref_crop_window(coords, lo, hi, delta, coarsen, buffer, align)
        got = crop_indices(coords, lo, hi, delta, coarsen, buffer, align)
        assert got == ref, (lo, hi, coarsen, buffer, align, got, ref)


@pytest.mark.python
def test_crop_indices_no_overlap_returns_none():
    coords = numpy.arange(0.0, 10.0, 0.5)
    assert crop_indices(coords, 20.0, 30.0, 0.5) is None   # request above range
    assert crop_indices(coords, -30.0, -20.0, 0.5) is None  # request below range


@pytest.mark.python
def test_crop_indices_subcell_window_raises():
    """A window inside the range but narrower than a cell is empty, not full.

    Distinguished from a genuine non-overlap (which returns None and has a
    documented full-file fall-back): quietly widening a crop the user narrowed
    too far is never right, so it raises.  EmptyCropWindow subclasses
    ValueError so existing handlers keep working.
    """
    coords = numpy.arange(0.0, 10.0, 0.5)
    assert issubclass(EmptyCropWindow, ValueError)
    with pytest.raises(EmptyCropWindow, match="between grid points"):
        crop_indices(coords, 1.1, 1.4, 0.5)   # strictly between 1.0 and 1.5


@pytest.mark.python
def test_axis_file_slice_ascending_is_the_identity_window():
    coords = numpy.arange(0.0, 10.0, 1.0)   # ascending
    sl, subset, flip = axis_file_slice(coords, False, 2, 8, 2, coords.size)
    assert flip is False
    assert sl == slice(2, 8, 2)
    numpy.testing.assert_array_equal(subset, coords[2:8:2])


@pytest.mark.python
def test_axis_file_slice_descending_selects_the_same_samples():
    """A descending file axis must yield the SAME physical samples.

    NetCDF/xarray lazy indexing requires a positive step, so a descending axis
    (latitude N->S) needs its ascending window remapped rather than reversed
    with a negative stride.  Pinned by comparing against the ascending view
    directly: same values, just in file order, with flip=True to say so.
    """
    asc = numpy.arange(0.0, 10.0, 1.0)
    desc = asc[::-1]                         # as stored in the file
    n = desc.size

    for lo, hi, step in [(0, n, 1), (2, 8, 2), (1, 9, 3), (3, 4, 1)]:
        sl, subset, flip = axis_file_slice(desc, True, lo, hi, step, n)
        assert flip is True
        assert sl.step > 0, "lazy indexing needs a positive stride"
        numpy.testing.assert_array_equal(subset, desc[sl])
        # Flipping the file-order subset must reproduce the ascending window.
        numpy.testing.assert_array_equal(subset[::-1], asc[lo:hi:step])


# ---------------------------------------------------------------------------
# coordinate_system wrap gate
# ---------------------------------------------------------------------------

@pytest.mark.python
def test_classify_lon_axis_positive_evidence():
    # Positive geographic evidence.
    assert classify_lon_axis(units="degrees_east") == "geographic"
    assert classify_lon_axis(units="degrees") == "geographic"
    assert classify_lon_axis(std_name="longitude") == "geographic"
    # Positive projected evidence.
    assert classify_lon_axis(units="m") == "projected"
    assert classify_lon_axis(units="km") == "projected"
    assert classify_lon_axis(std_name="projection_x_coordinate") == "projected"
    assert classify_lon_axis(lon_min=0.0, lon_max=500000.0) == "projected"
    # No metadata, in-range values -> ambiguous.
    assert classify_lon_axis(lon_min=0.0, lon_max=10.0) == "unknown"
    assert classify_lon_axis() == "unknown"


@pytest.mark.python
def test_resolve_wrap_geographic_run():
    # coordinate_system == 2 (lon-lat): wrapping allowed unless file is projected.
    assert resolve_wrap(2, "geographic") is True
    assert resolve_wrap(2, "unknown") is True
    with pytest.raises(ValueError, match="mismatch"):
        resolve_wrap(2, "projected")


@pytest.mark.python
def test_resolve_wrap_cartesian_run():
    # coordinate_system == 1 (Cartesian): never wrap; a positively geographic
    # file is a hard error, an ambiguous one is not.
    assert resolve_wrap(1, "projected") is False
    assert resolve_wrap(1, "unknown") is False
    with pytest.raises(ValueError, match="mismatch"):
        resolve_wrap(1, "geographic")


@pytest.mark.python
def test_resolve_wrap_unknown_system_is_legacy():
    # coordinate_system None -> per-file heuristic, no cross-check / no raise.
    assert resolve_wrap(None, "geographic") is True
    assert resolve_wrap(None, "unknown") is True
    assert resolve_wrap(None, "projected") is False


@pytest.mark.python
def test_is_geographic_lon_consistent_with_classify():
    # is_geographic_lon is 'not projected' -- the wrap-mechanics view.
    for units, sname, lo, hi in [
        (None, None, 0.0, 10.0),          # unknown -> geographic-for-wrap
        ("degrees_east", None, 0.0, 359.0),
        ("m", None, 0.0, 10.0),           # projected
        (None, None, 0.0, 500000.0),      # projected by range
    ]:
        expected = classify_lon_axis(units, sname, lo, hi) != "projected"
        assert is_geographic_lon(lo, hi, units, sname) is expected
