#!/usr/bin/env python
# encoding: utf-8
"""Pin the ``force_dry`` mask index convention at every site that reads it.

The mask is a 0/1 raster written by ``topotools.Topography.write`` with
``topo_type=3``.  Two facts about that file fix the convention the Fortran must
honour, and neither is documented anywhere else:

* the header's ``xlower``/``ylower`` are ``x[0]``, ``y[0]``, the **cell center**
  of the SW-most point, not a corner (``topotools.py`` ``grid_registration =
  'lower'``), and
* the data is written ``numpy.flipud(Z)``, so **file row 1 is the north row**,
  which ``read_force_dry`` stores verbatim as ``force_dry(:,1)``.

So for a patch whose cells align with the mask (the only case any of the lookups
run in), the correct map from a cell center to the stored array is::

    ii = floor((x - xlow_fdry) / dx_fdry) + 1
    jj = my_fdry - floor((y - ylow_fdry) / dy_fdry)

These tests assert the resulting dry/wet map directly, against an expectation
computed in Python from the same mask.  There is no stored golden: the oracle is
the convention above, restated independently in ``_expected_dry``.

Two failure modes are deliberately separated, because a fix can close one and
leave the other:

1. **The shift.**  Tests A/C/D use asymmetric, corner-touching masks, so a
   one-cell shift in either axis, or a transpose, is a hard failure.
2. **The halo.**  Test B uses a mask strictly smaller than the domain.  Fortran
   ``int()`` truncates toward zero rather than flooring, so a cell just outside
   the mask yields a quotient in ``(-1, 0)`` and lands on column 1 or row
   ``my_fdry`` instead of failing the ``(ii>=1) .and. ...`` bounds test.  Only
   ``floor()`` rejects it.  Test B is the one that catches that.

The ``src/2d/bouss/`` copies of these lookups cannot be exercised here (the
bouss build needs PETSc and MPI, which CI does not install).  Their parity with
the ``shallow`` copies is asserted in ``tests/test_force_dry_convention.py``.
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.test as gtest
import clawpack.geoclaw.topotools as topotools
from clawpack.pyclaw import solution

testdir = Path(__file__).parent

TOPO_Z = -10.0          # flat bathymetry
H0 = 10.0               # ambient depth, sea_level = 0

# ``data.py`` writes ForceDry.tend with "%.3f", so anything below 0.0005 s
# silently becomes 0.000 and disables the mask.  Never set this from tfinal.
TEND = 1.0e9


def _flat_topo(path):
    """Flat topotype-3 bathymetry covering (and overhanging) the domain."""
    topo = topotools.Topography(topo_func=lambda x, y: TOPO_Z + 0.0 * x)
    topo.topo_type = 3
    topo.x = np.linspace(-1.0, 11.0, 25)
    topo.y = np.linspace(-1.0, 9.0, 21)
    topo.write(path, topo_type=3, Z_format="%22.15e")


def _write_mask(path, x, y, Z):
    """Write a 0/1 force_dry mask in the documented format.

    ``set_xyZ`` requires ``Z.shape == (len(y), len(x))`` with ``y`` ascending and
    ``Z[0]`` the south row; ``write`` then flips it so the file is north-first.
    """
    mask = topotools.Topography()
    mask.set_xyZ(np.asarray(x), np.asarray(y), np.asarray(Z, dtype=float))
    mask.write(path, topo_type=3, Z_format="%1i")
    return np.asarray(Z)


def _expected_dry(state, x0, y0, delta, Z):
    """The dry/wet map the mask *should* produce on this patch.

    Independent restatement of the convention: a patch cell center sits exactly
    on a mask cell center (both lookups only run when the resolutions match), so
    the mask index is the nearest integer offset from ``(x0, y0)``, which are
    themselves cell centers.  Cells with no mask cell are untouched.
    """
    X, Y = state.grid.c_centers
    i = np.rint((X - x0) / delta).astype(int)
    j = np.rint((Y - y0) / delta).astype(int)
    ny, nx = Z.shape
    inside = (i >= 0) & (i < nx) & (j >= 0) & (j < ny)
    dry = np.zeros(X.shape, dtype=bool)
    dry[inside] = Z[j[inside], i[inside]] == 1
    return dry


@pytest.fixture(scope="module")
def plain_xgeoclaw(tmp_path_factory):
    """Build ``xgeoclaw`` once for the whole module.

    Every case here differs only in runtime data, so one ``make new`` replaces
    four rebuilds of ~140 Fortran sources.
    """
    build_dir = tmp_path_factory.mktemp("force_dry_build")
    builder = gtest.GeoClawTestRunner(build_dir, test_path=testdir)
    builder.build_executable()
    return build_dir / builder.executable_name


def _install_executable(runner, prebuilt):
    """Place a shared prebuilt executable where ``run_code`` expects it."""
    shutil.copy(prebuilt, runner.temp_path / runner.executable_name)


def _run(tmp_path, prebuilt, mask_x, mask_y, mask_Z, **setrun_kwargs):
    """Write inputs, run, and return (runner, mask origin, mask spacing, Z)."""
    topo_path = tmp_path / "flat.tt3"
    mask_path = tmp_path / "force_dry.tt3"
    _flat_topo(topo_path)
    Z = _write_mask(mask_path, mask_x, mask_y, mask_Z)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data(topo_path=str(topo_path), mask_path=str(mask_path),
                    tend_force_dry=TEND, **setrun_kwargs)
    runner.write_data()
    _install_executable(runner, prebuilt)
    runner.run_code()

    delta = float(mask_x[1] - mask_x[0])
    return runner, float(mask_x[0]), float(mask_y[0]), delta, Z


def _assert_exact_dry_map(state, x0, y0, delta, Z, context=""):
    """Frame 0 is written before any time stepping, so h is exactly 0 or H0."""
    expected = _expected_dry(state, x0, y0, delta, Z)
    h = np.asarray(state.q[0])
    actual = h == 0.0
    if not np.array_equal(actual, expected):
        X, Y = state.grid.c_centers
        wrong = np.argwhere(actual != expected)
        detail = ", ".join(
            f"(x={X[tuple(k)]:g}, y={Y[tuple(k)]:g}): "
            f"got {'dry' if actual[tuple(k)] else 'wet'}, "
            f"want {'dry' if expected[tuple(k)] else 'wet'}"
            for k in wrong[:8])
        pytest.fail(f"force_dry mask landed on the wrong cells{context}. "
                    f"{len(wrong)} cell(s) differ: {detail}")
    np.testing.assert_allclose(h[~expected], H0, rtol=1e-13)


@pytest.mark.regression
@pytest.mark.tsunami
def test_force_dry_qinit_asymmetric_mask(tmp_path, plain_xgeoclaw):
    """qinit.f90 at t0: an asymmetric mask pins the index map exactly.

    The mask covers the whole domain and is 1 in four cells chosen so that every
    way of getting the map wrong is a visible failure:

    * NW corner      -> the westmost column and the north-edge row flip,
    * SE corner      -> the eastmost column and the south row,
    * two interior cells with ``i != j`` -> a transposed map.

    On an unfixed tree the NW cell is never dried at all (``ii`` underflows to 0
    and is rejected by the bounds test) and the other three are displaced.
    """
    x = np.arange(10) + 0.5          # 0.5 .. 9.5, dx = 1
    y = np.arange(8) + 0.5           # 0.5 .. 7.5
    Z = np.zeros((8, 10))
    Z[7, 0] = 1                      # NW corner
    Z[0, 9] = 1                      # SE corner
    Z[5, 2] = 1                      # interior, i != j
    Z[2, 6] = 1                      # interior, i != j

    runner, x0, y0, delta, Z = _run(tmp_path, plain_xgeoclaw, x, y, Z,
                                    amr_levels_max=1, total_steps=1)

    sol = solution.Solution(0, path=runner.temp_path)
    _assert_exact_dry_map(sol.states[0], x0, y0, delta, Z)


@pytest.mark.regression
@pytest.mark.tsunami
def test_force_dry_outside_mask_is_untouched(tmp_path, plain_xgeoclaw):
    """qinit.f90: cells outside the mask rectangle must never be forced dry.

    This is the ``int()`` versus ``floor()`` test.  The mask is all ones and
    strictly interior, so the expected dry set is exactly the mask rectangle.
    Truncation toward zero pulls the column immediately west and the row
    immediately south of the mask into the mask's edge cells, which no bounds
    test can catch because the resulting indices are in range.
    """
    x = np.arange(4) + 3.5           # 3.5 .. 6.5
    y = np.arange(4) + 2.5           # 2.5 .. 5.5
    Z = np.ones((4, 4))

    runner, x0, y0, delta, Z = _run(tmp_path, plain_xgeoclaw, x, y, Z,
                                    amr_levels_max=1, total_steps=1)

    sol = solution.Solution(0, path=runner.temp_path)
    state = sol.states[0]

    expected = _expected_dry(state, x0, y0, delta, Z)
    assert expected.sum() == 16, "test setup: expected a 4x4 dry block"

    dry = np.asarray(state.q[0]) == 0.0
    leaked = dry & ~expected
    if leaked.any():
        X, Y = state.grid.c_centers
        where = ", ".join(f"(x={X[tuple(k)]:g}, y={Y[tuple(k)]:g})"
                          for k in np.argwhere(leaked)[:8])
        pytest.fail(
            f"{leaked.sum()} cell(s) outside the force_dry rectangle were "
            f"forced dry: {where}. Fortran int() truncates toward zero, so a "
            f"cell west or south of the mask maps into it; use floor().")

    _assert_exact_dry_map(state, x0, y0, delta, Z)


@pytest.mark.regression
@pytest.mark.tsunami
def test_force_dry_qinit_on_refined_patch(tmp_path, plain_xgeoclaw):
    """qinit.f90 on a level-2 patch, whose origin is not the domain origin.

    At t0 every level is initialized through ``qinit`` (``amrclaw`` ``ginit.f``),
    so a flagregion active from t1 = 0 gives a refined patch at t0.  Catches any
    error that depends on patch origin rather than domain origin.
    """
    x = np.arange(3.25, 6.0, 0.5)    # 3.25 .. 5.75, dx = 0.5 (level 2)
    y = np.arange(3.25, 6.0, 0.5)
    Z = np.zeros((len(y), len(x)))
    Z[0, 0] = 1                      # SW corner of the mask
    Z[-1, 1] = 1                     # north row, off-diagonal
    Z[2, -1] = 1                     # east column, off-diagonal

    runner, x0, y0, delta, Z = _run(
        tmp_path, plain_xgeoclaw, x, y, Z,
        amr_levels_max=2, total_steps=1,
        flagregion=[2, 2, 0.0, 1.0e9, 2.0, 7.0, 2.0, 7.0])

    sol = solution.Solution(0, path=runner.temp_path)
    fine = [s for s in sol.states if s.patch.level == 2]
    assert fine, "no level-2 patch was created; the flagregion did not fire"

    for state in fine:
        _assert_exact_dry_map(state, x0, y0, delta, Z,
                              context=" on a level-2 patch at t0")


@pytest.mark.regression
@pytest.mark.tsunami
def test_force_dry_filval_on_late_regrid(tmp_path, plain_xgeoclaw):
    """filval.f90: a patch created by a regrid, not by initialization.

    ``qinit`` handles every level at t0, so ``filval`` is only reached from
    ``regrid`` (``amrclaw`` ``gfixup.f``).  The flagregion therefore starts at
    ``t1 > 0``: no level 2 exists at t0, and the first regrid builds one through
    ``filval``.  The mask is written at level-2 resolution so the ``ddxy`` guard
    leaves level 1 alone and the coarse solution stays at rest.

    The run uses a fixed dt far below the CFL limit, so cells dried at the
    regrid refill by only a fraction of H0 over the whole run and the dry/wet
    map is still unambiguous several steps later.
    """
    x = np.arange(3.25, 6.0, 0.5)
    y = np.arange(3.25, 6.0, 0.5)
    Z = np.zeros((len(y), len(x)))
    Z[0, 0] = 1
    Z[-1, 1] = 1
    Z[2, -1] = 1

    dt = 1.0e-4
    total_steps = 8
    runner, x0, y0, delta, Z = _run(
        tmp_path, plain_xgeoclaw, x, y, Z,
        amr_levels_max=2, total_steps=total_steps, dt=dt,
        flagregion=[2, 2, 0.5 * dt, 1.0e9, 2.0, 7.0, 2.0, 7.0])

    sol = solution.Solution(total_steps, path=runner.temp_path)
    fine = [s for s in sol.states if s.patch.level == 2]
    assert fine, "no level-2 patch was created; the flagregion did not fire"

    for state in fine:
        expected = _expected_dry(state, x0, y0, delta, Z)
        h = np.asarray(state.q[0])
        actual = h < 0.5 * H0
        if not np.array_equal(actual, expected):
            X, Y = state.grid.c_centers
            wrong = np.argwhere(actual != expected)
            detail = ", ".join(
                f"(x={X[tuple(k)]:g}, y={Y[tuple(k)]:g}, h={h[tuple(k)]:g})"
                for k in wrong[:8])
            pytest.fail(
                f"force_dry mask landed on the wrong cells in filval "
                f"({len(wrong)} cell(s)): {detail}")
        assert (h[~expected] > 0.9 * H0).all(), \
            "cells outside the mask lost more depth than the run should allow"
