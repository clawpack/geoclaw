# encoding: utf-8
"""Minimal Cartesian setrun for the force_dry mask-index oracle.

The ``force_dry`` mask is applied by a pure index lookup from a cell center into
a raster: ``qinit.f90`` does it at ``t0``, ``filval.f90``/``filpatch.f90`` do it
for patches created later.  Nothing about that lookup depends on the physics, so
a tiny flat-bathymetry Cartesian domain is enough to pin the index convention
exactly, with no bathymetry gradient or wave to confound the dry/wet map.

Flat ``Z = -10`` with ``sea_level = 0`` gives an ambient ``h = 10`` everywhere,
so a cell the mask forces dry is exactly ``h = 0`` in frame 0 (written before
any time stepping).  The dry/wet map is then a clean boolean.

Parameterized by keyword:
- ``topo_path``       : absolute path to the flat topo file (topotype 3).
- ``mask_path``       : absolute path to the force_dry mask (topotype 3, 0/1).
- ``tend_force_dry``  : mask is applied while ``t <= tend``.  See the note in
  the test module: ``data.py`` writes this with ``"%.3f"``, so small values
  silently round to ``0.000`` and disable the mask.  Use a large value.
- ``amr_levels_max``  : 1 (qinit only) or 2 (level-2 patches).
- ``flagregion``      : optional ``[minlevel, maxlevel, t1, t2, x1, x2, y1, y2]``.
  A ``t1 > 0`` region is how the ``filval`` path is reached at all: at ``t0``
  every level is initialized through ``qinit``, so level 2 must not exist until
  the first regrid.
"""

from clawpack.clawutil import data


def setrun(claw_pkg="geoclaw", topo_path=None, mask_path=None,
           tend_force_dry=1.0e9, amr_levels_max=1, flagregion=None,
           dt=1.0e-4, total_steps=1):

    assert claw_pkg.lower() == "geoclaw", "Expected claw_pkg = 'geoclaw'"

    rundata = data.ClawRunData(claw_pkg, 2)
    clawdata = rundata.clawdata

    # --- Spatial domain: dx = dy = 1, so mask cells and cell centers align ---
    clawdata.num_dim = 2
    clawdata.lower[0] = 0.0
    clawdata.upper[0] = 10.0
    clawdata.lower[1] = 0.0
    clawdata.upper[1] = 8.0
    clawdata.num_cells[0] = 10
    clawdata.num_cells[1] = 8

    # --- System size: Cartesian, no friction, so topo is the only aux ---
    clawdata.num_eqn = 3
    clawdata.num_aux = 1
    clawdata.capa_index = 0

    # --- Time ---
    clawdata.t0 = 0.0
    clawdata.restart = False

    # Fixed-step output so a post-regrid frame lands at a predictable step.
    clawdata.output_style = 3
    clawdata.output_step_interval = 1
    clawdata.total_steps = total_steps
    clawdata.output_t0 = True
    clawdata.output_format = "ascii"
    clawdata.output_q_components = "all"
    clawdata.output_aux_components = "none"
    clawdata.verbosity = 0

    # --- Time stepping ---
    # Fixed dt, deliberately far below the CFL limit.  A cell dried at a regrid
    # then refills by only ~0.02 m per step against an ambient 10 m, so a frame
    # taken several steps after the regrid still carries an unambiguous signal.
    clawdata.dt_variable = False
    clawdata.dt_initial = dt
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

    # --- AMR ---
    amrdata = rundata.amrdata
    amrdata.amr_levels_max = amr_levels_max
    amrdata.refinement_ratios_x = [2]
    amrdata.refinement_ratios_y = [2]
    amrdata.refinement_ratios_t = [2]
    amrdata.aux_type = ["center"]
    amrdata.flag_richardson = False
    amrdata.flag2refine = (amr_levels_max > 1)
    amrdata.regrid_interval = 1
    amrdata.regrid_buffer_width = 2
    amrdata.verbosity_regrid = 0

    # Only the flagregion may flag; the solution is at rest otherwise.
    rundata.refinement_data.wave_tolerance = 1.0e9
    rundata.refinement_data.speed_tolerance = [1.0e9]
    rundata.refinement_data.variable_dt_refinement_ratios = False

    if flagregion is not None:
        from clawpack.amrclaw.data import FlagRegion
        region = FlagRegion(num_dim=2)
        region.convert_old_region(flagregion)
        rundata.flagregiondata.flagregions.append(region)

    # --- GeoClaw geometry / topo ---
    geo_data = rundata.geo_data
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 1      # Cartesian
    geo_data.coriolis_forcing = False
    geo_data.friction_forcing = False
    geo_data.dry_tolerance = 1.0e-3
    geo_data.sea_level = 0.0

    rundata.topo_data.topofiles.append([3, topo_path])

    # --- force_dry ---
    from clawpack.geoclaw.data import ForceDry
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.variable_eta_init = False
    force_dry = ForceDry()
    force_dry.tend = tend_force_dry
    force_dry.fname = mask_path
    rundata.qinit_data.force_dry_list.append(force_dry)

    return rundata


if __name__ == "__main__":
    import sys
    setrun(*sys.argv[1:]).write()
