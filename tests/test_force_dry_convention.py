#!/usr/bin/env python
# encoding: utf-8
"""The force_dry file contract, and parity across the five Fortran lookups.

Two guards that need no Fortran build:

1. ``test_force_dry_file_contract`` writes a mask exactly the way the docs and
   the ``eta_init_force_dry`` example do, and asserts the two properties the
   Fortran lookup depends on: the header origin is a **cell center**, and the
   data is written **north row first**.  Neither is stated anywhere in the
   Python source, the docs, or the file itself, which is why the Fortran drifted
   away from it.

2. ``test_force_dry_index_map_parity`` asserts the five sites that read the mask
   agree with each other.  They are hand-maintained copies of one expression;
   they disagreed for years.  The two ``src/2d/bouss/`` copies cannot be
   exercised by a test run at all (that build needs PETSc and MPI, which CI does
   not install), so this is their only coverage.
"""

import re
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.topotools as topotools

geoclaw_root = Path(__file__).parent.parent

# Every site that maps a cell center into the force_dry array.
FORCE_DRY_SITES = [
    Path("src/2d/shallow/qinit.f90"),
    Path("src/2d/shallow/filpatch.f90"),
    Path("src/2d/shallow/filval.f90"),
    Path("src/2d/bouss/filpatch.f90"),
    Path("src/2d/bouss/filval.f90"),
]


@pytest.mark.python
def test_force_dry_file_contract(tmp_path):
    """A written mask is cell-center registered and north-row first.

    The Fortran reads the third and fourth header lines as ``xlow_fdry`` and
    ``ylow_fdry`` and stores the data rows verbatim, so these two properties are
    precisely what makes ``ii = floor(...) + 1`` and ``jj = my_fdry - floor(...)``
    the right map.  If either ever changes, the Fortran must change with it.
    """
    x = np.arange(5) + 0.5            # 0.5 .. 4.5, dx = 1
    y = np.arange(3) + 0.5            # 0.5 .. 2.5
    Z = np.zeros((3, 5))
    Z[0, 0] = 1                       # SW corner
    Z[2, 4] = 1                       # NE corner
    Z[1, 3] = 1                       # asymmetric interior

    mask = topotools.Topography()
    mask.set_xyZ(x, y, Z)
    path = tmp_path / "force_dry.tt3"
    mask.write(path, topo_type=3, Z_format="%1i")

    lines = path.read_text().strip().split("\n")
    header, rows = lines[:6], lines[6:]

    def value(line):
        return float(line.split()[0])

    assert int(value(header[0])) == 5, "ncols"
    assert int(value(header[1])) == 3, "nrows"
    assert value(header[2]) == pytest.approx(x[0]), (
        "header origin must be the cell CENTER of the SW-most point (x[0]), "
        "not a corner; the Fortran +1 on the column index depends on it")
    assert value(header[3]) == pytest.approx(y[0]), "ylower is y[0]"
    assert value(header[4]) == pytest.approx(1.0), "cellsize"

    data = np.array([[int(v) for v in row.split()] for row in rows])
    assert data.shape == Z.shape
    np.testing.assert_array_equal(data[0], Z[-1], err_msg=(
        "first data row must be the NORTH row; read_force_dry stores rows "
        "verbatim, so force_dry(:,1) is north and the jj flip depends on it"))
    np.testing.assert_array_equal(data, np.flipud(Z))


@pytest.mark.python
def test_force_dry_index_map_parity():
    """All five lookups use the same expression, and none of them uses int().

    ``int()`` truncates toward zero, so a cell just west or south of the mask
    produces a quotient in ``(-1, 0)``, lands on an in-range index, and is
    silently forced dry.  ``floor()`` is what rejects it.
    """
    column = re.compile(
        r"ii\s*=\s*floor\(\s*\(\s*\w+\s*-\s*xlow_fdry\s*\+\s*1d-7\s*\)"
        r"\s*/\s*dx_fdry\s*\)\s*\+\s*1")
    row = re.compile(
        r"jj\s*=\s*my_fdry\s*-\s*floor\(\s*\(\s*\w+\s*-\s*ylow_fdry"
        r"\s*\+\s*1d-7\s*\)\s*/\s*dy_fdry\s*\)")

    problems = []
    for rel in FORCE_DRY_SITES:
        path = geoclaw_root / rel
        assert path.exists(), f"missing force_dry lookup site: {rel}"
        source = path.read_text()
        if not column.search(source):
            problems.append(f"{rel}: column index is not "
                            "'floor((x - xlow_fdry + 1d-7) / dx_fdry) + 1'")
        if not row.search(source):
            problems.append(f"{rel}: row index is not "
                            "'my_fdry - floor((y - ylow_fdry + 1d-7) / dy_fdry)'")
        for line in source.split("\n"):
            if "fdry" in line and re.search(r"\bint\s*\(", line):
                problems.append(f"{rel}: int() in a force_dry lookup: "
                                f"{line.strip()}")

    assert not problems, (
        "the force_dry index map differs across its copies:\n  "
        + "\n  ".join(problems))
