# encoding: utf-8
"""Minimal single-grid (no-AMR) setrun for the ``topo0save`` dtopo regression.

The case under test is *multiple topo files plus a dtopo*: ``topo0save``
decides which topo files have their t=0 state stashed in ``topo0work``, and
``topo_update`` skips any file whose ``topo0save`` is zero, so a file wrongly
marked zero never receives the seafloor deformation.

Everything here is deliberately inert -- flat bathymetry, still water, no
storm forcing, no AMR -- so the only thing that can move the bathymetry is the
dtopo.  The topo files are passed in as a list so one setrun serves both the
single-file control and the two-file case.
"""

from clawpack.clawutil import data
from clawpack.geoclaw import topotools


def setrun(claw_pkg="geoclaw", topo_paths=None, dtopo_path=None):

    assert claw_pkg.lower() == "geoclaw", "Expected claw_pkg = 'geoclaw'"

    rundata = data.ClawRunData(claw_pkg, 2)
    clawdata = rundata.clawdata

    # --- Spatial domain (small, fixed single grid) ---
    clawdata.num_dim = 2
    clawdata.lower[0] = -2.0
    clawdata.upper[0] = 2.0
    clawdata.lower[1] = -2.0
    clawdata.upper[1] = 2.0
    clawdata.num_cells[0] = 40
    clawdata.num_cells[1] = 40

    # --- System size ---
    clawdata.num_eqn = 3
    clawdata.num_aux = 3
    clawdata.capa_index = 2

    # --- Time ---
    clawdata.t0 = 0.0
    clawdata.restart = False

    # Frame 0 at t=0 (before the dtopo rise completes) and frame 1 well after
    # it, so the deformation is unambiguously present in the second frame.
    clawdata.output_style = 2
    clawdata.output_times = [0.0, 20.0]
    clawdata.output_format = "ascii"
    clawdata.output_q_components = "all"
    clawdata.output_aux_components = "all"   # need aux(1) = bathymetry
    clawdata.output_aux_onlyonce = False     # aux at every output time
    clawdata.verbosity = 0

    # --- Time stepping ---
    clawdata.dt_variable = True
    clawdata.dt_initial = 0.1
    clawdata.dt_max = 1e99
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0
    clawdata.steps_max = 100000

    # --- Method ---
    clawdata.order = 2
    clawdata.dimensional_split = "unsplit"
    clawdata.transverse_waves = 2
    clawdata.num_waves = 3
    clawdata.limiter = ["mc", "mc", "mc"]
    clawdata.use_fwaves = True
    clawdata.source_split = "godunov"
    clawdata.num_ghost = 2

    clawdata.bc_lower[0] = "extrap"
    clawdata.bc_upper[0] = "extrap"
    clawdata.bc_lower[1] = "extrap"
    clawdata.bc_upper[1] = "extrap"

    # --- AMR off: one fixed grid, so aux maps straight onto the domain ---
    amrdata = rundata.amrdata
    amrdata.amr_levels_max = 1
    amrdata.refinement_ratios_x = [2]
    amrdata.refinement_ratios_y = [2]
    amrdata.refinement_ratios_t = [2]
    amrdata.aux_type = ["center", "capacity", "yleft"]
    amrdata.flag_richardson = False
    amrdata.flag2refine = False
    amrdata.verbosity_regrid = 0

    # --- GeoClaw geometry / topo ---
    geo_data = rundata.geo_data
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2      # lat-lon
    geo_data.earth_radius = 6367.5e3
    geo_data.coriolis_forcing = False
    geo_data.friction_forcing = False
    geo_data.dry_tolerance = 1.0e-3
    geo_data.sea_level = 0.0

    for topo_path in topo_paths:
        topo = topotools.Topography()
        topo.path = str(topo_path)
        topo.topo_type = 3
        rundata.topo_data.topofiles.append(topo)

    rundata.dtopo_data.dtopofiles.append([3, str(dtopo_path)])
    rundata.dtopo_data.dt_max_dtopo = 1.0

    return rundata


if __name__ == "__main__":
    import sys
    setrun(*sys.argv[1:]).write()
