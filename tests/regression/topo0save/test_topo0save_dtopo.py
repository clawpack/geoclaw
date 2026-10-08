#!/usr/bin/env python
# encoding: utf-8
"""Regression: seafloor deformation must reach every topo file it overlaps.

``topo0save(i)`` marks the topo files whose t=0 state is stashed in
``topo0work``; ``topo_update`` (``src/2d/shallow/topo_update.f90``) skips any
file with ``topo0save <= 0``, so a file wrongly marked zero silently never
receives the dtopo deformation -- the seafloor does not move and no tsunami is
generated over that file's footprint.

``read_dtopo_settings`` sets the flag by intersecting each topo file against
the ``topo_for_dtopo`` entries, which occupy slots
``mtopofiles+1 ... mtopofiles+num_dtopo``.  With a **single** topo file a
wrong index range degenerates into comparing that file against itself, which
always "overlaps", so the flag comes out right for the wrong reason -- every
single-topo-file example passes regardless.  The bug only surfaces with two or
more topo files, which is why this case uses two.

The two topo files here are deliberately **disjoint** in longitude, mirroring
the antimeridian wrap-split that first exposed this: a global NetCDF DEM
cropped across +/-180 is written to ``topo.data`` as two separate entries.
They are separated by a gap narrow enough to stay inside the 1% domain
coverage tolerance in ``topo_module.f90`` (a wider gap is a hard ``stop``).

The dtopo is parameterized onto each half in turn, so the test does not depend
on which file happens to be listed first after the priority sort.
"""

from pathlib import Path

import numpy as np
import pytest

from clawpack.geoclaw import dtopotools, topotools
import clawpack.geoclaw.test as gtest
from clawpack.pyclaw import solution

testdir = Path(__file__).parent

# Flat seafloor everywhere, so any change in aux(1) is the dtopo and nothing
# else.
BASE_DEPTH = -1000.0
UPLIFT = 5.0

# The two topo files tile the domain apart from a narrow seam at x=0.  The gap
# must be non-zero (abutting files count as overlapping, which would mask the
# bug) but small: topo_module.f90 does `stop` when the uncovered area exceeds
# 1% of the domain.  0.01 deg across a 4 deg domain is 0.25%.
SEAM = 0.005
DOMAIN = (-2.0, 2.0, -2.0, 2.0)

# The topo files must be *finer* than the dtopo.  GeoClaw builds an internal
# "topo_for_dtopo" grid matching the dtopo and inserts it into the priority
# order by cell area; that grid is always deformed, so if it were the finest
# thing over the footprint it would supply the deformed seafloor no matter what
# topo0save said about the real topo files, and the bug would be invisible.
# Making the real topo finer puts it on top, so whether *it* got deformed is
# what the solution actually sees.
TOPO_DX = 0.02
DTOPO_DX = 0.1

# dtopo footprints, one well inside each half, clear of the seam and the
# domain edges.
FOOTPRINTS = {
    "west": (-1.3, -0.7),
    "east": (0.7, 1.3),
}


def _write_half_topos(tmp_path):
    """Two flat topotype-3 files, disjoint in x, together covering the domain."""
    paths = []
    for name, (x0, x1) in (("west", (-2.5, -SEAM)), ("east", (SEAM, 2.5))):
        topo = topotools.Topography(
            topo_func=lambda x, y: BASE_DEPTH + 0.0 * x)
        topo.topo_type = 3
        topo.x = np.arange(x0, x1 + 0.5 * TOPO_DX, TOPO_DX)
        topo.y = np.arange(-2.5, 2.5 + 0.5 * TOPO_DX, TOPO_DX)
        path = tmp_path / f"{name}.tt3"
        topo.write(path, topo_type=3, Z_format="%22.15e")
        paths.append(path)
    return paths


def _write_full_topo(tmp_path):
    """Single flat topotype-3 file covering the whole domain (the control)."""
    topo = topotools.Topography(topo_func=lambda x, y: BASE_DEPTH + 0.0 * x)
    topo.topo_type = 3
    topo.x = np.arange(-2.5, 2.5 + 0.5 * TOPO_DX, TOPO_DX)
    topo.y = np.arange(-2.5, 2.5 + 0.5 * TOPO_DX, TOPO_DX)
    path = tmp_path / "full.tt3"
    topo.write(path, topo_type=3, Z_format="%22.15e")
    return [path]


def _write_dtopo(tmp_path, half):
    """A uniform +UPLIFT step over a small box inside the requested half.

    Single-time (``mt = 1``) so the displacement is instantaneous: GeoClaw's
    ``topo_update`` takes its ``mtdtopo == 1`` branch, applies the full ``dZ``
    once with no time interpolation, and immediately finalizes.  A dtopo that
    ramps over an interval instead couples the answer to the time-stepping
    (the deformation gets frozen at whatever fraction the first step reached),
    which has nothing to do with what this test is checking.
    """
    x0, x1 = FOOTPRINTS[half]
    x = np.arange(x0, x1 + 0.5 * DTOPO_DX, DTOPO_DX)
    y = np.arange(-0.3, 0.3 + 0.5 * DTOPO_DX, DTOPO_DX)
    X, Y = np.meshgrid(x, y)

    dtopo = dtopotools.DTopography()
    dtopo.x, dtopo.y, dtopo.X, dtopo.Y = x, y, X, Y
    dtopo.times = [0.0]
    dtopo.dZ = np.array([UPLIFT * np.ones_like(X)])

    path = tmp_path / f"uplift_{half}.tt3"
    dtopo.write(path, dtopo_type=3)
    return path


def _run(tmp_path, topo_paths, dtopo_path, prebuilt):
    """Write inputs, install the shared build, run, return the final aux."""
    import shutil

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data(topo_paths=[str(p) for p in topo_paths],
                    dtopo_path=str(dtopo_path))
    runner.write_data()
    shutil.copy(prebuilt, runner.temp_path / runner.executable_name)
    runner.run_code()

    # Frame 1 is t=20, well after the dtopo finishes rising.
    sol = solution.Solution(1, path=runner.temp_path, read_aux=True)
    state = sol.states[0]
    return state, np.asarray(state.aux[0])


# Sample well inside the dtopo footprint.  The deformation is uniform there,
# so the expected bathymetry is exact; sampling right up to the footprint edge
# instead picks up cells where the dtopo has been bilinearly interpolated onto
# the finer topo grid and averaged into domain cells, which smears the value by
# a few percent and would force a tolerance loose enough to be uninformative.
INSET = 0.15


def _bathymetry_over(state, bathy, half):
    """Mean bathymetry over the interior of the dtopo footprint for *half*."""
    x0, x1 = FOOTPRINTS[half]
    xc, yc = state.grid.c_centers
    mask = ((xc >= x0 + INSET) & (xc <= x1 - INSET)
            & (yc >= -0.3 + INSET) & (yc <= 0.3 - INSET))
    assert mask.any(), "footprint mask selected no cells"
    return bathy[mask].mean()


@pytest.fixture(scope="module")
def xgeoclaw(tmp_path_factory):
    """Build ``xgeoclaw`` once for the whole module."""
    build_dir = tmp_path_factory.mktemp("topo0save_build")
    builder = gtest.GeoClawTestRunner(build_dir, test_path=testdir)
    builder.build_executable()
    return build_dir / builder.executable_name


@pytest.mark.regression
@pytest.mark.parametrize("half", ["west", "east"])
def test_dtopo_reaches_both_topo_files(tmp_path, xgeoclaw, half):
    """Deformation must be applied whichever of two topo files it lands in.

    Before the fix the second-listed file was compared against itself (always
    "overlapping") while the first was compared against it (disjoint, so
    marked as needing no update), and the ``west`` case saw no deformation at
    all -- the seafloor stayed at its initial depth.
    """
    topo_paths = _write_half_topos(tmp_path)
    dtopo_path = _write_dtopo(tmp_path, half)

    state, bathy = _run(tmp_path, topo_paths, dtopo_path, xgeoclaw)
    got = _bathymetry_over(state, bathy, half)

    assert got == pytest.approx(BASE_DEPTH + UPLIFT, abs=1e-6), (
        f"dtopo over the {half} topo file was not applied: bathymetry is "
        f"{got}, expected {BASE_DEPTH + UPLIFT}. topo0save almost certainly "
        f"marked that file as needing no update."
    )


@pytest.mark.regression
@pytest.mark.parametrize("half", ["west", "east"])
def test_two_topo_files_match_single_file_control(tmp_path, xgeoclaw, half):
    """Splitting the topography into two files must not change the answer.

    The single-file control exercises the degenerate self-comparison that
    always worked; requiring the two-file case to match it pins the two
    against each other rather than against a hand-computed number.
    """
    dtopo_path = _write_dtopo(tmp_path, half)

    ctl_dir = tmp_path / "control"
    ctl_dir.mkdir()
    ctl_state, ctl_bathy = _run(
        ctl_dir, _write_full_topo(ctl_dir), dtopo_path, xgeoclaw)

    split_dir = tmp_path / "split"
    split_dir.mkdir()
    split_state, split_bathy = _run(
        split_dir, _write_half_topos(split_dir), dtopo_path, xgeoclaw)

    assert _bathymetry_over(split_state, split_bathy, half) == pytest.approx(
        _bathymetry_over(ctl_state, ctl_bathy, half), abs=1e-6)
