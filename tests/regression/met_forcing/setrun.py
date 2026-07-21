# encoding: utf-8
"""Minimal single-grid (no-AMR) met-forcing setrun for the Phase 0 aux oracle.

The storm forcing wind (u, v) and pressure aux fields are a pure function of
(x, y, t) via ``set_storm_fields`` — independent of the evolving solution — so a
tiny flat-bathymetry domain with AMR disabled is sufficient to characterize the
forcing contract exactly, with no AMR interpolation or solver feedback to
confound a diff.

Parameterized by keyword so one build serves both forcing families:
- ``forcing``   : ``"holland80"`` (parametric) or ``"data"`` (gridded netCDF).
- ``topo_path`` : absolute path to a flat topo file (topotype 3).
- ``storm_path``: absolute path to the storm file / gridded descriptor.
"""

from clawpack.clawutil import data


def setrun(claw_pkg="geoclaw", forcing="holland80", topo_path=None,
           storm_path=None, amr_levels_max=1, wind_refine=False,
           R_refine=False):

    assert claw_pkg.lower() == "geoclaw", "Expected claw_pkg = 'geoclaw'"

    rundata = data.ClawRunData(claw_pkg, 2)
    clawdata = rundata.clawdata

    # --- Spatial domain (small, fixed single grid) ---
    clawdata.num_dim = 2
    clawdata.lower[0] = -5.0
    clawdata.upper[0] = 5.0
    clawdata.lower[1] = 15.0
    clawdata.upper[1] = 25.0
    clawdata.num_cells[0] = 20
    clawdata.num_cells[1] = 20

    # --- System size ---
    clawdata.num_eqn = 3
    # 3 GeoClaw geometry slots + 1 friction + 3 storm (wind_u, wind_v, pressure)
    clawdata.num_aux = 3 + 1 + 3
    clawdata.capa_index = 2

    # --- Time ---
    clawdata.t0 = 0.0
    clawdata.restart = False

    # A few fixed output times covering the storm window.
    clawdata.output_style = 2
    clawdata.output_times = [0.0, 3.0 * 3600.0, 6.0 * 3600.0]
    clawdata.output_format = "ascii"       # ascii -> deterministic aux dump
    clawdata.output_q_components = "all"
    clawdata.output_aux_components = "all"   # dump fort.aNNNN
    clawdata.output_aux_onlyonce = False     # aux at every output time
    clawdata.verbosity = 0

    # --- Time stepping ---
    clawdata.dt_variable = True
    clawdata.dt_initial = 1.0
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

    # Boundary conditions: extrapolation on all sides.
    clawdata.bc_lower[0] = "extrap"
    clawdata.bc_upper[0] = "extrap"
    clawdata.bc_lower[1] = "extrap"
    clawdata.bc_upper[1] = "extrap"

    # --- AMR (single fixed grid by default; F1 test drives amr_levels_max=2) ---
    amrdata = rundata.amrdata
    amrdata.amr_levels_max = amr_levels_max
    amrdata.refinement_ratios_x = [2]
    amrdata.refinement_ratios_y = [2]
    amrdata.refinement_ratios_t = [2]
    amrdata.aux_type = ["center", "capacity", "yleft",
                        "center", "center", "center", "center"]
    amrdata.flag_richardson = False
    # Use flag2refine2 (wind/R_refine criteria) when refining.
    amrdata.flag2refine = (amr_levels_max > 1)
    amrdata.verbosity_regrid = 0

    # --- GeoClaw geometry / topo ---
    geo_data = rundata.geo_data
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2      # lat-lon
    geo_data.earth_radius = 6367.5e3
    geo_data.coriolis_forcing = True
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.dry_tolerance = 1.0e-2
    geo_data.sea_level = 0.0

    rundata.topo_data.topofiles.append([3, topo_path])

    # --- Storm / met forcing ---
    surge_data = rundata.surge_data
    surge_data.wind_forcing = True
    surge_data.drag_law = 1
    surge_data.pressure_forcing = True
    surge_data.wind_index = 4          # 0-based -> Fortran 5 (wind_u), 6 (wind_v)
    surge_data.pressure_index = 6      # 0-based -> Fortran 7 (pressure)
    surge_data.display_landfall_time = False
    surge_data.wind_refine = wind_refine
    surge_data.R_refine = R_refine
    if forcing == "data":
        surge_data.storm_specification_type = "data"
    else:
        surge_data.storm_specification_type = forcing   # e.g. "holland80"
    surge_data.storm_file = storm_path

    return rundata


if __name__ == "__main__":
    import sys
    setrun(*sys.argv[1:]).write()
