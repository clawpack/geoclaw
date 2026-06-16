"""Pytest regression tests for Fortran-side topography preprocessing.

These exercise the per-file preprocessing attributes (``z_shift``,
``negate_z``, ``x_shift``, ``crop_extent``, ``buffer``, ``coarsen``,
``align``) end-to-end through
the real GeoClaw executable: each variant writes topo file(s) whose stored
data has the *inverse* transformation baked in, so the preprocessing applied
by Fortran at read time must reconstruct the same physical bowl bathymetry
as the unmodified baseline.  All variants therefore compare against a single
shared reference gauge series.

Every variant is constructed so that a silently broken preprocessing step
produces *physically different* bathymetry near the gauge, not just a no-op:

- ``z_shift`` / ``negate_z`` / ``x_shift``: the file data is wrong by the
  inverse transformation, so skipping the step leaves wrong bathymetry
  everywhere.
- ``z_shift`` also plants a no-data cell at the gauge location whose
  ``topo_missing`` replacement equals the true bowl value; if the Fortran
  fill-cell guard fails, the shift digs a 5 m hole under the gauge.
- ``crop`` / ``crop_buffer``: the cropped file is *corrupted* east of the
  crop boundary and a second file supplies the correct data there.  If the
  crop (or its buffer arithmetic) is not applied, the corrupted or uncovered
  region reaches the gauge.
- ``coarsen`` / ``coarsen_align``: the fine file is corrupted at every point
  the subsampling should skip, so a missed or mis-indexed stride keeps
  corrupted points.

The initial condition (``qinit.f90``) sets ``h = max(0, eta - b)`` directly
from the topo aux array, so bathymetry errors are immediately visible in the
gauge depth series.

To (re)generate the reference data, run the baseline variant with ``--save``::

    pytest test_bowl_slosh_preprocessing.py -k baseline --save

Scenarios from ``NOTES_topotools_preprocessing.md`` (Prompt 5):
z_shift (incl. fill-cell guard), negate_z (+ double-negate identity),
priority order, and the ``topo_type=1`` preprocessing guard; plus coverage
of the Fortran ``crop_extent``/``buffer``/``x_shift`` read paths.
"""

from pathlib import Path
import subprocess

import numpy as np
import pytest

import clawpack.geoclaw.test as test
import clawpack.geoclaw.topotools as topotools

# Bowl geometry (matches qinit.f90)
A = 1.0
H0 = 0.1

GAUGE_ID = 11
GAUGE_XY = (0.5, 0.5)

# Grid spacing shared by all topo files so that grid points always align
# with the baseline grid (crop snapping is then exact).
DX = 0.02

# Matches Topography.no_data_value default, written to the file header.
NO_DATA = -99999.0


def _bowl(x, y):
    return H0 * (x**2 + y**2) / A**2 - H0


def _corrupt_east_of(seam):
    """Bowl west of (and at) the seam; 0.05 too deep east of it.

    The threshold sits half a cell east of the seam so that the seam grid
    point itself is never corrupted by floating-point round-off (the seam
    line survives the crop and must hold correct values).
    """
    return lambda x, y: _bowl(x, y) - 0.05 * (x > seam + DX / 2)


def _write_topo_file(path, zfunc=_bowl, x=None, y=None, topo_type=2):
    """Write a topo file; defaults produce the exact baseline bowl."""
    topo = topotools.Topography(topo_func=zfunc)
    topo.x = np.linspace(-2.0, 2.0, 201) if x is None else x
    topo.y = np.linspace(-2.0, 2.0, 201) if y is None else y
    topo.write(path, topo_type=topo_type, Z_format="%22.15e")


def _write_fine_corrupt_file(path, x, y):
    """Write a fine-resolution bowl corrupted off the even-even lattice.

    Points whose row or column index is odd are 0.05 too deep; subsampling
    with coarsen=2 from an even start index keeps only correct points, while
    any indexing error lands on corrupted ones.
    """
    X, Y = np.meshgrid(x, y)
    Z = _bowl(X, Y)
    corrupt = np.ones(Z.shape, dtype=bool)
    corrupt[::2, ::2] = False
    Z[corrupt] -= 0.05
    topo = topotools.Topography()
    topo.x = x
    topo.y = y
    topo.Z = Z
    topo.write(path, topo_type=3, Z_format="%22.15e")


def _topo_entry(path, topo_type=2, **attrs):
    """Construct a Topography topofiles entry with preprocessing attributes."""
    t = topotools.Topography()
    t.path = path
    t.topo_type = topo_type
    for name, value in attrs.items():
        setattr(t, name, value)
    return t


def _setup_baseline(temp_path):
    _write_topo_file(temp_path / "bowl.topotype2")
    return [_topo_entry("bowl.topotype2")], None


def _setup_z_shift(temp_path):
    # File stores elevations offset by +5; z_shift removes it at read time.
    # The grid point at the gauge location is marked no-data: Fortran replaces
    # it with topo_missing, which we set to the TRUE bowl value there, and the
    # fill-cell guard must then skip the z_shift for that cell.  If the guard
    # fails, the gauge sits over a 5 m hole and the comparison fails.
    x = np.linspace(-2.0, 2.0, 201)
    y = np.linspace(-2.0, 2.0, 201)
    X, Y = np.meshgrid(x, y)
    Z = _bowl(X, Y) + 5.0
    i_gauge = 125                       # x[125] = y[125] = 0.5, the gauge point
    Z[i_gauge, i_gauge] = NO_DATA
    topo = topotools.Topography()
    topo.x = x
    topo.y = y
    topo.Z = Z
    topo.write(temp_path / "bowl.topotype2", topo_type=2, Z_format="%22.15e")

    return ([_topo_entry("bowl.topotype2", z_shift=-5.0)],
            _bowl(*GAUGE_XY))


def _setup_negate_z(temp_path):
    # File stores depths (positive down); negate_z restores elevations
    # independently of the legacy negative-topo_type convention.
    _write_topo_file(temp_path / "bowl.topotype2",
                     zfunc=lambda x, y: -_bowl(x, y))
    return [_topo_entry("bowl.topotype2", negate_z=True)], None


def _setup_double_negate(temp_path):
    # Negative topo_type (legacy flip) plus negate_z: the two sign flips
    # cancel, so the file stores the true elevations.
    _write_topo_file(temp_path / "bowl.topotype2")
    return [_topo_entry("bowl.topotype2", topo_type=-2, negate_z=True)], None


def _setup_x_shift(temp_path):
    # File grid is registered 0.5 too far west; x_shift slides it back.
    _write_topo_file(temp_path / "bowl.topotype2",
                     zfunc=lambda x, y: _bowl(x + 0.5, y),
                     x=np.linspace(-2.5, 1.5, 201))
    return [_topo_entry("bowl.topotype2", x_shift=0.5)], None


def _setup_crop(temp_path):
    # The first-listed (oversized) file is corrupted east of x=0 and cropped
    # to the western half; a second file supplies the correct eastern half.
    # If the crop is not applied, the corrupted region covers the gauge.
    _write_topo_file(temp_path / "west.topotype2",
                     zfunc=_corrupt_east_of(0.0),
                     x=np.linspace(-3.0, 3.0, 301),
                     y=np.linspace(-3.0, 3.0, 301))
    _write_topo_file(temp_path / "east.topotype2",
                     x=np.linspace(0.0, 2.0, 101))
    return ([_topo_entry("west.topotype2", crop_extent=[-2.0, 0.0, -2.0, 2.0]),
             _topo_entry("east.topotype2")],
            None)


def _setup_crop_buffer(temp_path):
    # Same idea with the seam at x=0.26, two cells (buffer=2) inside the
    # crop_extent on every side: the buffer must grow the cropped region to
    # exactly [-2, 0.26] x [-2, 2] to meet the eastern file.  If the buffer
    # is ignored, an uncovered sliver opens 0.24 west of the gauge (within
    # wave-travel distance for t_final = 0.5).
    seam = 0.26
    _write_topo_file(temp_path / "west.topotype2",
                     zfunc=_corrupt_east_of(seam),
                     x=np.linspace(-3.0, 3.0, 301),
                     y=np.linspace(-3.0, 3.0, 301))
    _write_topo_file(temp_path / "east.topotype2",
                     x=np.linspace(seam, 2.0, 88))
    extent = [-2.0 + 2 * DX, seam - 2 * DX, -2.0 + 2 * DX, 2.0 - 2 * DX]
    return ([_topo_entry("west.topotype2", crop_extent=extent, buffer=2),
             _topo_entry("east.topotype2")],
            None)


def _setup_coarsen(temp_path):
    # File at twice the baseline resolution, corrupted off the even-even
    # lattice; coarsen=2 keeps exactly the baseline grid points.  If the
    # subsampling is skipped or mis-indexed, corrupted points reach the
    # gauge.  (Written as topo_type=3 to exercise the row-format path.)
    _write_fine_corrupt_file(temp_path / "fine.topotype3",
                             x=np.linspace(-2.0, 2.0, 401),
                             y=np.linspace(-2.0, 2.0, 401))
    return [_topo_entry("fine.topotype3", topo_type=3, coarsen=2)], None


def _setup_coarsen_align(temp_path):
    # Fine grid registered so that the naive subsample start lands on a
    # corrupted (odd-index) point: the crop lower edge (-2.015) selects
    # x = -2.01 (index 1) as the first in-window point, and align=(0.5, 0.5)
    # must shift the start to x = -2.00 (index 2, even lattice), recovering
    # exactly the baseline grid.  If align is ignored, the kept lattice is
    # entirely corrupted.
    _write_fine_corrupt_file(temp_path / "fine.topotype3",
                             x=np.linspace(-2.02, 2.0, 403),
                             y=np.linspace(-2.02, 2.0, 403))
    return ([_topo_entry("fine.topotype3", topo_type=3, coarsen=2,
                         crop_extent=[-2.015, 2.0, -2.015, 2.0],
                         align=(0.5, 0.5))],
            None)


def _setup_priority(temp_path):
    # Two same-resolution files covering the same extent.  The last-listed
    # file wins in the overlap (reversed mtopoorder: last = highest priority;
    # the Python priority sort is stable for equal cell areas), so listing the
    # true bowl last makes the run match the baseline.
    _write_topo_file(temp_path / "bowl.topotype2")
    _write_topo_file(temp_path / "deep.topotype2",
                     zfunc=lambda x, y: _bowl(x, y) - 0.05)
    return ([_topo_entry("deep.topotype2"), _topo_entry("bowl.topotype2")],
            None)


_VARIANTS = {
    "baseline": _setup_baseline,
    "z_shift": _setup_z_shift,
    "negate_z": _setup_negate_z,
    "double_negate": _setup_double_negate,
    "x_shift": _setup_x_shift,
    "crop": _setup_crop,
    "crop_buffer": _setup_crop_buffer,
    "coarsen": _setup_coarsen,
    "coarsen_align": _setup_coarsen_align,
    "priority": _setup_priority,
}


def _run_bowl(tmp_path, topofiles, topo_missing=None):
    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    runner.set_data()

    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.5

    runner.rundata.topo_data.topofiles = topofiles
    if topo_missing is not None:
        runner.rundata.topo_data.topo_missing = topo_missing

    runner.rundata.gaugedata.gauges = []
    runner.rundata.gaugedata.gauges.append(
        [GAUGE_ID, GAUGE_XY[0], GAUGE_XY[1], 0.0, 1e10])

    runner.write_data()
    runner.build_executable()
    runner.run_code()
    return runner


@pytest.mark.regression
@pytest.mark.parametrize("variant", list(_VARIANTS))
def test_preprocessing_identity(tmp_path: Path, save: bool, variant: str):
    """Each preprocessing variant must reproduce the baseline gauge series."""
    topofiles, topo_missing = _VARIANTS[variant](tmp_path)
    runner = _run_bowl(tmp_path, topofiles, topo_missing=topo_missing)

    # Only the baseline writes reference data, so --save cannot replace the
    # reference with another variant's (tolerance-equal, but not bitwise
    # identical) output.
    runner.check_gauge(gauge_id=GAUGE_ID, indices=(2, 3),
                       rtol=1e-4, atol=1e-4,
                       save=(save and variant == "baseline"))


@pytest.mark.regression
def test_priority_order_detectable(tmp_path: Path):
    """Negative control for the priority variant.

    With the deeper file listed last it must win in the overlap, producing
    gauge output that *differs* from the baseline reference.  This guards
    against the priority test passing vacuously (e.g. if overlapping files
    were averaged or the first file always won).
    """
    _write_topo_file(tmp_path / "bowl.topotype2")
    _write_topo_file(tmp_path / "deep.topotype2",
                     zfunc=lambda x, y: _bowl(x, y) - 0.05)
    topofiles = [_topo_entry("bowl.topotype2"), _topo_entry("deep.topotype2")]
    runner = _run_bowl(tmp_path, topofiles)

    with pytest.raises(AssertionError):
        runner.check_gauge(gauge_id=GAUGE_ID, indices=(2, 3),
                           rtol=1e-4, atol=1e-4)


@pytest.mark.regression
def test_topo_type1_preprocessing_guard(tmp_path: Path):
    """topo_type=1 plus any preprocessing attribute must stop the Fortran run.

    The guard in read_topo_settings prints a diagnostic and executes
    ``stop 1``, which surfaces here as a CalledProcessError from runclaw.
    """
    _write_topo_file(tmp_path / "bowl.topotype1", topo_type=1)
    topofiles = [_topo_entry("bowl.topotype1", topo_type=1, z_shift=5.0)]

    with pytest.raises(subprocess.CalledProcessError):
        _run_bowl(tmp_path, topofiles)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
