#!/usr/bin/env python
# encoding: utf-8
"""End-to-end regression for the topo descriptor-crop path (``topo_type=4``).

This is the first test that compiles GeoClaw and feeds it a NetCDF topo file
through a descriptor written by ``TopoInspector.topo_entries()``.  That path had
no coverage at all, which is how the bug below shipped.

``read_topo_file`` selects the file subset three ways, in priority order:
descriptor ``crop_bounds``, then ``topo_crop_extent``, then the AMR domain.  The
buffer expansion (``nbuf4``) was applied on the second branch only, so any file
carrying ``crop_bounds`` -- which is *every* file ``topo_entries()`` writes --
silently discarded its declared ``buffer``.

The failure is not subtle once you look for it: GeoClaw refuses to start with
"topo arrays do not cover domain" for a domain the user correctly sized against
crop-plus-buffer.  The two cases here are the before/after of exactly that.

Each case asserts on the topo grid GeoClaw reports in ``fort.geo`` -- the
window actually loaded -- rather than only on whether the run survived, so an
over-buffered or ignored crop fails as loudly as a dropped one.

NetCDF is required; the whole module skips without it.
"""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.test as gtest
import clawpack.geoclaw.topotools as topotools
from clawpack.geoclaw import netcdf_utils as ncutils
from clawpack.geoclaw.data import TopographyData

testdir = Path(__file__).parent

# Topo file, generously larger than any crop below.  DELTA is a negative power
# of two so every coordinate is exact in binary: with 0.05 the file's own
# spacing came out as 4.9999999999997e-2 and the crop_bounds comparison
# (an inclusive >= / <= on the raw coordinates) landed one index differently on
# each edge, which made the window asymmetric for reasons having nothing to do
# with what is under test.
TOPO_X = (-100.0, -80.0)
TOPO_Y = (20.0, 30.0)
DELTA = 0.125

# The crop window written into the descriptor, and the buffer requested.  All
# four edges are integer multiples of DELTA away from the file origin.
CROP = (-95.0, -90.0, 22.0, 25.0)
BUFFER = 4
BUFFER_DEG = BUFFER * DELTA          # 0.5 degrees on each side

# Margin used when placing a test domain relative to a topo window.  The domain
# must not merely *touch* the topo edge: GeoClaw's coverage test is an
# inequality on floats, and a domain whose bounds coincide exactly with the topo
# bounds passes it and then crashes intermittently in the boundary
# interpolation.  Two cells of slack keeps these tests about the buffer.
MARGIN = 2 * DELTA

pytestmark = pytest.mark.regression


def _netcdf_build_vars():
    """make_vars enabling a NetCDF build, or None if unavailable.

    Same probe as the met_forcing suite: ``nf-config`` for the Fortran flags,
    plus ``nc-config --libs`` because ``nf-config --flibs`` names ``-lnetcdf``
    without its ``-L`` path.
    """
    nf = shutil.which("nf-config")
    nc = shutil.which("nc-config")
    if nf is None:
        return None
    try:
        fflags = subprocess.check_output([nf, "--fflags"], text=True).strip()
        flibs = subprocess.check_output([nf, "--flibs"], text=True).strip()
        if nc is not None:
            flibs += " " + subprocess.check_output(
                [nc, "--libs"], text=True).strip()
    except (subprocess.CalledProcessError, OSError):
        return None
    return {"USE_NETCDF": "1", "NETCDF_FFLAGS": fflags,
            "NETCDF_LFLAGS": flibs}


@pytest.fixture(scope="module")
def netcdf_xgeoclaw(tmp_path_factory):
    """Build the NetCDF-enabled ``xgeoclaw`` once for the module."""
    make_vars = _netcdf_build_vars()
    if make_vars is None:
        pytest.skip("NetCDF (nf-config/nc-config) unavailable; the topo "
                    "descriptor path is topo_type=4 only.")
    build_dir = tmp_path_factory.mktemp("topo_crop_build")
    builder = gtest.GeoClawTestRunner(build_dir, test_path=testdir)
    builder.build_executable(make_vars=make_vars)
    return build_dir / builder.executable_name


def _write_nc_topo(path):
    """A CF-compliant NetCDF bathymetry file spanning TOPO_X x TOPO_Y.

    The field is a plain south-to-north ramp: nothing here tests values, only
    which subset of them Fortran loads.
    """
    netCDF4 = pytest.importorskip("netCDF4")

    x = np.arange(TOPO_X[0], TOPO_X[1] + 1e-9, DELTA)
    y = np.arange(TOPO_Y[0], TOPO_Y[1] + 1e-9, DELTA)
    Z = -100.0 + np.outer(np.linspace(0.0, 90.0, y.size), np.ones_like(x))

    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('lon', x.size)
        ds.createDimension('lat', y.size)
        v = ds.createVariable('lon', 'f8', ('lon',))
        v[:] = x
        v.units = 'degrees_east'
        v.standard_name = 'longitude'
        v = ds.createVariable('lat', 'f8', ('lat',))
        v[:] = y
        v.units = 'degrees_north'
        v.standard_name = 'latitude'
        v = ds.createVariable('elevation', 'f8', ('lat', 'lon'))
        v[:] = Z
        v.units = 'm'
        v.standard_name = 'height_above_mean_sea_level'
        v.positive = 'up'
        ds.Conventions = 'CF-1.8'
    return path


def _write_topo_data(out_path, nc_path, buffer):
    """Write topo.data via topo_entries(), so the descriptor has crop_bounds.

    Going through ``topo_entries()`` rather than setting ``crop_extent`` alone
    is the point: it is what puts a ``crop_bounds`` line in the descriptor, and
    therefore what selects the Fortran branch under test.
    """
    with ncutils.TopoInspector(str(nc_path), crop_bounds=CROP) as insp:
        entries = insp.topo_entries()

    topos = []
    for _topo_type, _path, meta in entries:
        t = topotools.Topography()
        t.path = str(nc_path)
        t.topo_type = 4
        t.crop_extent = list(CROP)
        t.buffer = buffer
        t._netcdf_meta = meta
        topos.append(t)

    td = TopographyData()
    td.topofiles = topos
    td.write(out_file=str(out_path))
    return len(entries)


def _ring_domain():
    """A domain lying in the buffer ring: outside CROP, inside CROP+buffer.

    Both gaps are MARGIN wide, so the case is decided by whether the buffer was
    applied and not by float comparisons at a coincident edge.
    """
    return (CROP[0] - MARGIN, CROP[1] + MARGIN,
            CROP[2] - MARGIN, CROP[3] + MARGIN)


def _topo_window(tmp_path):
    """Parse the topo grid GeoClaw actually loaded out of ``fort.geo``.

    read_topo_header writes ``mx = N  x = (lo,hi)`` / ``my = ...`` for each
    file; with one topo file the first pair is ours.
    """
    text = (tmp_path / "fort.geo").read_text()
    m = re.search(r"mx\s*=\s*(\d+)\s+x\s*=\s*\(\s*([-\d.E+]+)\s*,"
                  r"\s*([-\d.E+]+)\s*\)", text)
    n = re.search(r"my\s*=\s*(\d+)\s+y\s*=\s*\(\s*([-\d.E+]+)\s*,"
                  r"\s*([-\d.E+]+)\s*\)", text)
    assert m is not None and n is not None, f"no topo grid in fort.geo:\n{text}"
    return (float(m.group(2)), float(m.group(3)),
            float(n.group(2)), float(n.group(3)),
            int(m.group(1)), int(n.group(1)))


def _run(tmp_path, prebuilt, domain, buffer):
    """Set up and run one case; return (returncode, stdout)."""
    nc_path = tmp_path / "topo.nc"
    _write_nc_topo(nc_path)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data()
    cd = runner.rundata.clawdata
    cd.lower = [domain[0], domain[2]]
    cd.upper = [domain[1], domain[3]]
    # Cell count is immaterial to the topo-coverage check; keep it small but
    # non-degenerate.
    cd.num_cells = [20, 12]
    runner.write_data()

    # write_data() emits its own topo.data; overwrite it with the
    # descriptor version, which is what this test is actually about.
    n_entries = _write_topo_data(tmp_path / "topo.data", nc_path, buffer)
    assert n_entries == 1, "crop lies wholly inside the file; expected 1 entry"

    shutil.copy(prebuilt, tmp_path / runner.executable_name)
    (tmp_path / "_output").mkdir(exist_ok=True)
    proc = subprocess.run([str(tmp_path / runner.executable_name)],
                          cwd=tmp_path, capture_output=True, text=True)
    # Kept for post-mortem: pytest truncates a long assertion message.
    (tmp_path / "run.log").write_text(proc.stdout + proc.stderr)
    return proc.returncode, proc.stdout + proc.stderr


def test_descriptor_crop_honors_buffer(tmp_path, netcdf_xgeoclaw):
    """A domain inside crop+buffer but outside crop alone must run.

    This is the regression.  With ``nbuf4`` pinned to 0 on the ``crop_bounds``
    branch the loaded topo covered only CROP, leaving the outer 0.20-degree ring
    empty, and GeoClaw stopped with "topo arrays do not cover domain" -- an
    abort, for inputs that were correctly specified.
    """
    domain = _ring_domain()
    rc, out = _run(tmp_path, netcdf_xgeoclaw, domain, BUFFER)

    assert "topo arrays do not cover domain" not in out, (
        "buffer was not applied to the descriptor crop; the topo loaded "
        f"covers only {CROP}.\n{out}")
    assert rc == 0, out

    # Stronger than "it ran": the loaded window must be CROP grown by exactly
    # BUFFER points on every side.  A test that only checks for the abort would
    # also pass if the reader over-buffered, or ignored crop_bounds and loaded
    # the whole file.
    xll, xhi, yll, yhi, mx, my = _topo_window(tmp_path)
    assert xll == pytest.approx(CROP[0] - BUFFER_DEG)
    assert xhi == pytest.approx(CROP[1] + BUFFER_DEG)
    assert yll == pytest.approx(CROP[2] - BUFFER_DEG)
    assert yhi == pytest.approx(CROP[3] + BUFFER_DEG)
    assert mx == round((CROP[1] - CROP[0]) / DELTA) + 2 * BUFFER + 1
    assert my == round((CROP[3] - CROP[2]) / DELTA) + 2 * BUFFER + 1


def test_descriptor_crop_without_buffer_does_not_cover_ring(
        tmp_path, netcdf_xgeoclaw):
    """The complement: with buffer=0 the same domain is genuinely uncovered.

    Without this, the test above would still pass if the reader started
    ignoring ``crop_bounds`` entirely and loaded the whole file -- which would
    be a different bug with the same symptom.  Here the coverage failure is the
    correct answer, and seeing it proves the crop is being applied at all.
    """
    domain = _ring_domain()
    rc, out = _run(tmp_path, netcdf_xgeoclaw, domain, 0)

    assert "topo arrays do not cover domain" in out, (
        "buffer=0 should leave the ring outside crop_bounds uncovered; the "
        f"reader appears not to be applying crop_bounds at all.\n{out}")


def test_descriptor_crop_inside_window_runs_either_way(
        tmp_path, netcdf_xgeoclaw):
    """Control: a domain wholly inside CROP runs with no buffer at all.

    Pins that the descriptor crop itself is loaded correctly, independently of
    the buffer question.
    """
    inset = 4 * DELTA
    domain = (CROP[0] + inset, CROP[1] - inset,
              CROP[2] + inset, CROP[3] - inset)

    rc, out = _run(tmp_path, netcdf_xgeoclaw, domain, 0)

    assert "topo arrays do not cover domain" not in out, out
    assert rc == 0, out

    # The window is CROP itself: unbuffered, and not the whole file.
    xll, xhi, yll, yhi, _mx, _my = _topo_window(tmp_path)
    assert (xll, xhi, yll, yhi) == pytest.approx(
        (CROP[0], CROP[1], CROP[2], CROP[3]))


# ===========================================================================
# Antimeridian: two entries for one file, with different lon_wrap_offset
#
# This is the end-to-end wrap.  Nothing had ever fed GeoClaw two descriptor
# entries pointing at the same file with different lon_wrap_offset values --
# the Python side was tested, and the Fortran side was assumed.
# ===========================================================================

# A global file so a crop can genuinely run off its longitude extent.
WRAP_TOPO_X = (-180.0, 180.0)
WRAP_TOPO_Y = (-40.0, 0.0)
WRAP_DELTA = 0.25

# Continuous-coordinate crop spanning the seam: covered by 170..180 (offset
# -360) plus -180..-150 (offset 0).
WRAP_CROP = (-190.0, -150.0, -30.0, -10.0)


def _write_global_nc(path):
    """A CF-compliant global NetCDF topo file spanning the antimeridian.

    Depth varies with longitude so the two sides of the seam carry visibly
    different values; a constant field would hide a mis-stitched join.
    """
    netCDF4 = pytest.importorskip("netCDF4")

    x = np.arange(WRAP_TOPO_X[0], WRAP_TOPO_X[1] + 1e-9, WRAP_DELTA)
    y = np.arange(WRAP_TOPO_Y[0], WRAP_TOPO_Y[1] + 1e-9, WRAP_DELTA)
    Z = -3000.0 + 10.0 * np.cos(np.radians(x))[None, :] * np.ones((y.size, 1))

    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('lon', x.size)
        ds.createDimension('lat', y.size)
        v = ds.createVariable('lon', 'f8', ('lon',))
        v[:] = x
        v.units = 'degrees_east'
        v.standard_name = 'longitude'
        v = ds.createVariable('lat', 'f8', ('lat',))
        v[:] = y
        v.units = 'degrees_north'
        v.standard_name = 'latitude'
        v = ds.createVariable('elevation', 'f8', ('lat', 'lon'))
        v[:] = Z
        v.units = 'm'
        v.standard_name = 'height_above_mean_sea_level'
        v.positive = 'up'
        ds.Conventions = 'CF-1.8'
    return path


def _run_wrapped(tmp_path, prebuilt, domain, buffer=0):
    """Set up and run a case whose topo.data has two wrapped entries."""
    nc_path = tmp_path / "global.nc"
    _write_global_nc(nc_path)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data()
    cd = runner.rundata.clawdata
    cd.lower = [domain[0], domain[2]]
    cd.upper = [domain[1], domain[3]]
    cd.num_cells = [20, 12]
    runner.write_data()

    # Written the way a user would: one Topography with a continuous
    # cross-seam crop.  The writer expands it into two descriptor entries.
    topo = topotools.Topography()
    topo.path = str(nc_path)
    topo.topo_type = 4
    topo.crop_extent = list(WRAP_CROP)
    topo.buffer = buffer

    td = TopographyData()
    td.topofiles = [topo]
    td.write(out_file=str(tmp_path / "topo.data"))

    text = (tmp_path / "topo.data").read_text()
    assert text.count("lon_wrap_offset") == 2, (
        f"expected two wrapped entries, got:\n{text}")

    shutil.copy(prebuilt, tmp_path / runner.executable_name)
    (tmp_path / "_output").mkdir(exist_ok=True)
    proc = subprocess.run([str(tmp_path / runner.executable_name)],
                          cwd=tmp_path, capture_output=True, text=True)
    (tmp_path / "run.log").write_text(proc.stdout + proc.stderr)
    return proc.returncode, proc.stdout + proc.stderr


@pytest.mark.netcdf
def test_wrapped_entries_cover_a_cross_seam_domain(tmp_path, netcdf_xgeoclaw):
    """A domain straddling the antimeridian is covered by the two entries.

    The domain sits at -185..-155 in continuous coordinates, i.e. it spans the
    date line.  Neither entry covers it alone: the offset-0 entry supplies
    -180..-150 and the offset -360 entry supplies -190..-180 by reading the
    file's 170..180 strip.  If Fortran ignored lon_wrap_offset, or applied it
    before selecting on crop_bounds, the west half would be missing and the
    coverage check would fail.
    """
    domain = (-185.0, -155.0, -28.0, -12.0)
    rc, out = _run_wrapped(tmp_path, netcdf_xgeoclaw, domain)

    assert "topo arrays do not cover domain" not in out, (
        f"the wrapped pair did not cover a cross-seam domain.\n{out}")
    assert rc == 0, out


@pytest.mark.netcdf
def test_wrapped_entries_report_shifted_extents(tmp_path, netcdf_xgeoclaw):
    """Both entries must be reported in *domain* coordinates.

    fort.geo records each topo grid after lon_wrap_offset is applied, so the
    wrapped entry has to appear west of -180 rather than at its file position
    of +170.  Seeing +170 here would mean the offset never reached the grid.
    """
    domain = (-185.0, -155.0, -28.0, -12.0)
    _run_wrapped(tmp_path, netcdf_xgeoclaw, domain)

    text = (tmp_path / "fort.geo").read_text()
    windows = re.findall(
        r"mx\s*=\s*\d+\s+x\s*=\s*\(\s*([-\d.E+]+)\s*,\s*([-\d.E+]+)\s*\)",
        text)
    assert len(windows) == 2, f"expected two topo grids in fort.geo:\n{text}"

    lows = sorted(float(lo) for lo, _hi in windows)
    highs = sorted(float(hi) for _lo, hi in windows)
    # Wrapped entry shifted to -190..-180; unwrapped entry at -180..-150.
    assert lows[0] == pytest.approx(-190.0)
    assert highs[-1] == pytest.approx(-150.0)
    # Nothing left sitting at its unshifted file position.
    assert all(float(hi) < 0.0 for _lo, hi in windows), (
        f"an entry kept its file longitude instead of the domain one: "
        f"{windows}")
