#!/usr/bin/env python
# encoding: utf-8
r"""Phase 0 "nothing broke" canary for the meteorological-forcing refactor.

bowl-slosh is a pure tsunami case that uses no storm forcing, yet it compiles
the shallow ``src2`` / ``setaux`` / ``b4step2`` / ``flag2refine2`` sources that
``use storm_module``.  Its output must therefore be **byte-identical across
every storm-refactor commit** — the moment this canary moves, shared
(non-storm) code has been perturbed.

This is stricter than the tolerance-based gauge/fgmax checks in
``test_bowl_slosh.py``: it byte-compares the full output surface (``fort.q*``,
``fort.t*``, gauge files, and fgmax output) against committed goldens under
``regression_data/canary/``.  ASCII output makes the comparison portable on a
given platform/compiler.

Regenerate the goldens after an *intentional* change with ``GEOCLAW_REGEN=1``.
Byte-identity assumes a fixed platform/compiler; a diff means either shared code
changed (a real canary trip) or the build environment changed.
"""

import os
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.test as test
import clawpack.geoclaw.fgmax_tools as fgmax_tools
import clawpack.geoclaw.topotools as topotools

# Deterministic output only: solution frames, gauge series, and fgmax results.
# Timing/log/parameter files (timing.*, pyclaw.log, fort.debug) are intentionally
# excluded — they carry run-to-run timing or path noise.
CANARY_GLOBS = ("fort.q*", "fort.t*", "gauge*.txt", "fgmax*.txt")


def _golden_name(produced_name):
    """Map a produced output filename to its committed golden filename.

    ``fort.q*`` / ``fort.t*`` are excluded by common global gitignores
    (``fort.q*`` etc.), so store their goldens under a ``golden_`` prefix that
    escapes those patterns while staying obviously paired.  Other outputs
    (gauge*.txt, fgmax*.txt) are not ignored and keep their names.
    """
    if produced_name.startswith("fort."):
        return "golden_" + produced_name[len("fort."):]
    return produced_name


def _assert_output_bytes(temp_path, golden_dir):
    """Byte-compare the canary output surface against committed goldens."""
    temp_path = Path(temp_path)
    golden_dir = Path(golden_dir)
    regen = bool(os.environ.get("GEOCLAW_REGEN"))

    produced = sorted(
        {p for pat in CANARY_GLOBS for p in temp_path.glob(pat)},
        key=lambda p: p.name,
    )
    assert produced, "canary produced no output files to compare"

    if regen or not golden_dir.exists():
        golden_dir.mkdir(parents=True, exist_ok=True)
        for p in produced:
            (golden_dir / _golden_name(p.name)).write_bytes(p.read_bytes())
        if not regen:
            pytest.skip(f"Canary baseline created in {golden_dir.name} "
                        "(rerun to assert)")
        return

    mismatches = []
    expected_golden_names = set()
    for p in produced:
        gname = _golden_name(p.name)
        expected_golden_names.add(gname)
        golden = golden_dir / gname
        if not golden.exists():
            mismatches.append(f"{p.name}: no golden committed ({gname})")
        elif golden.read_bytes() != p.read_bytes():
            mismatches.append(f"{p.name}: bytes differ")
    # Also flag goldens that were expected but not produced.
    for golden in golden_dir.iterdir():
        if golden.name not in expected_golden_names:
            mismatches.append(f"{golden.name}: golden not reproduced")
    assert not mismatches, (
        "Canary output changed — shared (non-storm) code was perturbed:\n  "
        + "\n  ".join(mismatches)
        + "\nIf intentional, regenerate with GEOCLAW_REGEN=1."
    )


@pytest.mark.regression
@pytest.mark.tsunami
def test_bowl_slosh_canary(tmp_path: Path):
    """Non-surge bowl-slosh output must be byte-identical across the refactor."""
    runner = test.GeoClawTestRunner(tmp_path, test_path=Path(__file__).parent)
    runner.set_data()

    # Short, deterministic run.  ASCII output for portable byte comparison.
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.tfinal = 0.5
    runner.rundata.clawdata.output_format = "ascii"

    runner.rundata.gaugedata.gauges = [[1, 0.5, 0.5, 0, 1e10]]

    # fgmax grid (same as test_bowl_slosh) so fgmax output is exercised.
    runner.rundata.fgmax_data.num_fgmax_val = 2
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2
    fg.x1, fg.x2, fg.y1, fg.y2 = -2.0, 2.0, -2.0, 2.0
    fg.dx = 0.1
    fg.tstart_max = 0.0
    fg.tend_max = 1.0e10
    fg.dt_check = 0.1
    fg.min_level_check = 2
    fg.arrival_tol = 1.0e-2
    runner.rundata.fgmax_data.fgmax_grids.append(fg)

    runner.write_data()

    # Build the bowl-slosh topography (same as test_bowl_slosh).
    h0 = 0.1
    topo = topotools.Topography(
        topo_func=lambda x, y: h0 * (x**2 + y**2) - h0)
    topo.topo_type = 2
    topo.x = np.linspace(-2.0, 2.0, 200)
    topo.y = np.linspace(-2.0, 2.0, 200)
    topo.write(Path(runner.temp_path) / "bowl.topotype2", topo_type=2,
               Z_format="%22.15e")

    runner.build_executable()
    runner.run_code()

    _assert_output_bytes(runner.temp_path,
                         Path(__file__).parent / "regression_data" / "canary")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
