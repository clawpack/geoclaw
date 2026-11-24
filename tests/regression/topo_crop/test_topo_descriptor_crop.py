#!/usr/bin/env python
# encoding: utf-8
"""End-to-end regression for the topo descriptor-crop path (``topo_type=4``).

This is the first test that compiles GeoClaw and feeds it a NetCDF topo file
through a descriptor written by ``TopoInspector.topo_entries()``.  That path had
no coverage at all, which is how the bug below shipped.

``read_topo_file`` used to select the file subset three ways, in priority
order: descriptor ``crop_bounds``, then ``topo_crop_extent``, then the AMR
domain.  The buffer expansion (``nbuf4``) was applied on the second branch
only, so any file carrying ``crop_bounds`` -- which was *every* file
``topo_entries()`` wrote -- silently discarded its declared ``buffer``.
topo.data no longer has ``crop_bounds`` (dev/design/topo_input_roadmap.md
§2.8): the crop is ``topo_crop_extent`` for every entry, and the cases below
still pin that the buffer is honoured.

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
# spacing came out as 4.9999999999997e-2 and the crop comparison
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
    """Write topo.data from hand-built ``topo_entries()`` entries.

    The entries carry the crop only as their metadata's file-coordinate
    ``crop_bounds``; ``crop_extent`` is deliberately left unset.  topo.data has
    no ``crop_bounds`` line, so ``TopographyData.write`` must turn it into the
    entry's ``crop_extent`` -- otherwise the crop would be lost and Fortran
    would fall back to the AMR domain.
    """
    with ncutils.TopoInspector(str(nc_path), crop_bounds=CROP) as insp:
        entries = insp.topo_entries()

    topos = []
    for _topo_type, _path, meta in entries:
        t = topotools.Topography()
        t.path = str(nc_path)
        t.topo_type = 4
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

    This is the regression.  With ``nbuf4`` pinned to 0 on the old ``crop_bounds``
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
    # also pass if the reader over-buffered, or ignored the crop and loaded
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
    ignoring the crop entirely and loaded the whole file -- which would
    be a different bug with the same symptom.  Here the coverage failure is the
    correct answer, and seeing it proves the crop is being applied at all.
    """
    domain = _ring_domain()
    rc, out = _run(tmp_path, netcdf_xgeoclaw, domain, 0)

    assert "topo arrays do not cover domain" in out, (
        "buffer=0 should leave the ring outside the crop uncovered; the "
        f"reader appears not to be applying the crop at all.\n{out}")
    # The coverage failure is fatal and must say so in the exit status.
    assert rc != 0, out


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
    before selecting the crop, the west half would be missing and the
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


# ===========================================================================
# coarsen + align with no aligned point inside the crop (#761)
#
# The align step moves the window start onto the coarse lattice through
# `align`.  When the crop is narrower than `coarsen` points, that can move the
# start past the end of the window.  apply_align_buffer_coarsen used to carry
# on, keeping one point outside the crop (possibly outside the file); it must
# stop and name the file instead, for both reader paths.
# ===========================================================================

# Two grid points (-95.0, -94.875) in x.  coarsen=3 aligned to -94.75 keeps
# ..., -95.125, -94.75, ..., so neither is on the lattice.  y is wide and its
# align (22.0) lies on the grid, so only x is at fault.
EMPTY_ALIGN_CROP = (-95.0, -94.875, 22.0, 25.0)
EMPTY_ALIGN = (-94.75, 22.0)
EMPTY_ALIGN_COARSEN = 3


def _run_empty_align(tmp_path, prebuilt, topo_type):
    """Run with a crop whose aligned coarse lattice misses it entirely."""
    nc_path = _write_nc_topo(tmp_path / "topo.nc")
    if topo_type == 4:
        topo_path = nc_path
    else:
        src = topotools.Topography()
        src.read(str(nc_path), topo_type=4)
        topo_path = tmp_path / "topo.tt3"
        src.write(str(topo_path), topo_type=3)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data()
    cd = runner.rundata.clawdata
    # The stop comes while reading topo, before the coverage check, so the
    # domain only needs to be valid.
    cd.lower = [-96.0, 22.0]
    cd.upper = [-90.0, 25.0]
    cd.num_cells = [20, 12]
    runner.write_data()

    topo = topotools.Topography()
    topo.path = str(topo_path)
    topo.topo_type = topo_type
    topo.crop_extent = list(EMPTY_ALIGN_CROP)
    topo.coarsen = EMPTY_ALIGN_COARSEN
    topo.align = list(EMPTY_ALIGN)

    td = TopographyData()
    td.topofiles = [topo]
    td.write(out_file=str(tmp_path / "topo.data"))

    shutil.copy(prebuilt, tmp_path / runner.executable_name)
    (tmp_path / "_output").mkdir(exist_ok=True)
    proc = subprocess.run([str(tmp_path / runner.executable_name)],
                          cwd=tmp_path, capture_output=True, text=True)
    (tmp_path / "run.log").write_text(proc.stdout + proc.stderr)
    return proc.returncode, proc.stdout + proc.stderr, topo_path


@pytest.mark.parametrize("topo_type", [3, 4])
def test_crop_with_no_aligned_point_stops(tmp_path, netcdf_xgeoclaw,
                                          topo_type):
    """The reader stops, naming the file and axis, instead of inventing a
    point outside the crop.  Python's crop_indices raises for the same crop
    (EmptyAlignedCropWindow); this pins the Fortran mirror to agree."""
    rc, out, topo_path = _run_empty_align(tmp_path, netcdf_xgeoclaw,
                                          topo_type)

    assert rc != 0, out
    assert "crop contains no point of the coarsened grid" in out, out
    assert topo_path.name in out, out
    assert re.search(r"axis\s*=\s*x", out), out


# ===========================================================================
# Fortran backstop: a crop that misses its file (#760)
#
# TopographyData.write rejects a crop the file does not cover, so a missing
# crop only reaches Fortran through a hand-edited or stale topo.data.  Fortran
# must still stop and name the file.  Types 2/3 already did; type 4 turned the
# entry into a 0x0 grid with extents (0, 0) and kept it in the topo list.
# ===========================================================================

# Wholly east of the file (TOPO_X = -100..-80).
MISSING_CROP = (-70.0, -65.0, 22.0, 25.0)


def _run_missing_crop(tmp_path, prebuilt, topo_type):
    """Write a valid topo.data, then move the crop off the file by editing it,
    as a stale or hand-edited topo.data would."""
    nc_path = _write_nc_topo(tmp_path / "topo.nc")
    if topo_type == 4:
        topo_path = nc_path
    else:
        src = topotools.Topography()
        src.read(str(nc_path), topo_type=4)
        topo_path = tmp_path / "topo.tt3"
        src.write(str(topo_path), topo_type=3)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data()
    cd = runner.rundata.clawdata
    cd.lower = [-96.0, 22.0]
    cd.upper = [-90.0, 25.0]
    cd.num_cells = [20, 12]
    runner.write_data()

    topo = topotools.Topography()
    topo.path = str(topo_path)
    topo.topo_type = topo_type
    topo.crop_extent = list(CROP)
    td = TopographyData()
    td.topofiles = [topo]
    td.write(out_file=str(tmp_path / "topo.data"))

    missing = " ".join(repr(v) for v in MISSING_CROP)
    text = (tmp_path / "topo.data").read_text()
    text, n = re.subn(r"^.*# crop_extent.*$", f"{missing}   # crop_extent",
                      text, flags=re.M)
    assert n == 1, text
    (tmp_path / "topo.data").write_text(text)

    shutil.copy(prebuilt, tmp_path / runner.executable_name)
    (tmp_path / "_output").mkdir(exist_ok=True)
    proc = subprocess.run([str(tmp_path / runner.executable_name)],
                          cwd=tmp_path, capture_output=True, text=True)
    (tmp_path / "run.log").write_text(proc.stdout + proc.stderr)
    return proc.returncode, proc.stdout + proc.stderr, topo_path


@pytest.mark.parametrize("topo_type", [3, 4])
def test_crop_missing_the_file_stops(tmp_path, netcdf_xgeoclaw, topo_type):
    """A crop wholly outside its file stops the run, naming the file."""
    rc, out, topo_path = _run_missing_crop(tmp_path, netcdf_xgeoclaw,
                                           topo_type)

    assert rc != 0, out
    assert "does not overlap topo file" in out, out
    assert topo_path.name in out, out


# ===========================================================================
# crop_extent is in domain coordinates under x_shift/y_shift (§2.8)
#
# The same Topography settings must load the same region for ASCII and NetCDF.
# Types 2/3 subtract the shift before cropping; type 4 used to crop the file's
# own coordinates and shift afterwards, so it loaded CROP + shift and a domain
# inside CROP was reported as not covered.
#
# Checked on values, not on fort.geo: for types 2/3 fort.geo reports the file
# header before the crop, so it cannot show which window was used.  The field
# is linear, which GeoClaw's bilinear cell integral reproduces exactly, so each
# cell's B is the field at the cell centre, shifted back into file coordinates.
# ===========================================================================

SHIFT = (2.0, 1.0)      # whole multiples of DELTA


def _linear_field(x_file, y_file):
    return -100.0 + 2.0 * (x_file - TOPO_X[0]) + 3.0 * (y_file - TOPO_Y[0])


def _write_linear_topo(tmp_path, topo_type):
    x = np.arange(TOPO_X[0], TOPO_X[1] + 1e-9, DELTA)
    y = np.arange(TOPO_Y[0], TOPO_Y[1] + 1e-9, DELTA)
    X, Y = np.meshgrid(x, y)
    src = topotools.Topography()
    src.x, src.y, src.Z = x, y, _linear_field(X, Y)
    if topo_type == 4:
        netCDF4 = pytest.importorskip("netCDF4")
        path = tmp_path / "topo.nc"
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
            v[:] = src.Z
            v.units = 'm'
            ds.Conventions = 'CF-1.8'
        return path
    path = tmp_path / "topo.tt3"
    src.write(str(path), topo_type=3)
    return path


@pytest.mark.parametrize("topo_type", [3, 4])
def test_shifted_crop_is_in_domain_coordinates(tmp_path, netcdf_xgeoclaw,
                                               topo_type):
    """With a shift, the loaded topo is CROP in domain coordinates."""
    from clawpack.pyclaw import Solution

    topo_path = _write_linear_topo(tmp_path, topo_type)

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data()
    cd = runner.rundata.clawdata
    inset = 4 * DELTA                 # domain strictly inside CROP
    cd.lower = [CROP[0] + inset, CROP[2] + inset]
    cd.upper = [CROP[1] - inset, CROP[3] - inset]
    cd.num_cells = [16, 8]
    cd.output_t0 = True               # aux is written with frame 0
    cd.output_aux_components = 'all'
    cd.output_aux_onlyonce = True
    runner.write_data()

    topo = topotools.Topography()
    topo.path = str(topo_path)
    topo.topo_type = topo_type
    topo.crop_extent = list(CROP)
    topo.x_shift, topo.y_shift = SHIFT
    td = TopographyData()
    td.topofiles = [topo]
    td.write(out_file=str(tmp_path / "topo.data"))

    shutil.copy(netcdf_xgeoclaw, tmp_path / runner.executable_name)
    (tmp_path / "_output").mkdir(exist_ok=True)
    proc = subprocess.run([str(tmp_path / runner.executable_name)],
                          cwd=tmp_path, capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    (tmp_path / "run.log").write_text(out)

    assert "topo arrays do not cover domain" not in out, (
        f"type {topo_type} cropped in file coordinates: the topo loaded is "
        f"CROP shifted by {SHIFT}, not CROP.\n{out}")
    assert proc.returncode == 0, out

    sol = Solution(0, path=str(tmp_path), file_format='ascii', read_aux=True)
    state = sol.states[0]
    xc, yc = state.grid.p_centers
    dy = state.grid.delta[1]
    B = state.aux[0]
    np.testing.assert_allclose(
        B, _linear_field(xc - SHIFT[0], _cos_weighted_mean(yc, dy) - SHIFT[1]),
        rtol=0, atol=1e-8)


def _cos_weighted_mean(yc, dy):
    """Mean latitude of each cell weighted by cos(latitude), in degrees.

    On the sphere (coordinate_system 2) the topo cell average is weighted by
    area, i.e. by cos(latitude) (bilinearintegral.f), so a field linear in y
    averages to its value at this latitude rather than at the cell centre.
    """
    k = np.pi / 180.0
    lo, hi = yc - dy / 2, yc + dy / 2
    # Integral of y cos(k y) is y sin(k y) / k + cos(k y) / k**2.
    num = (hi * np.sin(k * hi) + np.cos(k * hi) / k
           - lo * np.sin(k * lo) - np.cos(k * lo) / k)
    den = np.sin(k * hi) - np.sin(k * lo)
    return num / den


# ===========================================================================
# An old topo.data with a descriptor crop_bounds line stops (§2.8)
#
# Older GeoClaw wrote crop_bounds and Fortran preferred it to crop_extent.  A
# hand-built entry could carry its crop in crop_bounds alone, so ignoring the
# line would load the wrong window without a message; the reader stops instead.
# ===========================================================================

def test_old_topo_data_with_crop_bounds_stops(tmp_path, netcdf_xgeoclaw):
    nc_path = _write_nc_topo(tmp_path / "topo.nc")

    runner = gtest.GeoClawTestRunner(tmp_path, test_path=testdir)
    runner.set_data()
    cd = runner.rundata.clawdata
    cd.lower = [-95.0, 22.0]
    cd.upper = [-90.0, 25.0]
    cd.num_cells = [20, 12]
    runner.write_data()

    topo = topotools.Topography()
    topo.path = str(nc_path)
    topo.topo_type = 4
    topo.crop_extent = list(CROP)
    td = TopographyData()
    td.topofiles = [topo]
    td.write(out_file=str(tmp_path / "topo.data"))

    # Put back the line an older writer emitted, just before the blank line
    # that ends the descriptor.
    old = " ".join(repr(v) for v in CROP)
    text = (tmp_path / "topo.data").read_text()
    text, n = re.subn(r"^(fill_action\s*=.*)$",
                      rf"\1\ncrop_bounds    = {old}", text, flags=re.M)
    assert n == 1, text
    (tmp_path / "topo.data").write_text(text)

    shutil.copy(netcdf_xgeoclaw, tmp_path / runner.executable_name)
    (tmp_path / "_output").mkdir(exist_ok=True)
    proc = subprocess.run([str(tmp_path / runner.executable_name)],
                          cwd=tmp_path, capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    (tmp_path / "run.log").write_text(out)

    assert proc.returncode != 0, out
    assert "Regenerate topo.data" in out, out
    assert nc_path.name in out, out
