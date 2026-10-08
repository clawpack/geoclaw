#!/usr/bin/env python
# encoding: utf-8
"""Run configuration for the topo descriptor-crop regression case.

The domain is set by the test, not here: each case picks bounds that sit in a
specific relationship to the topo file's crop window (inside it, or inside the
crop-plus-buffer ring), which is the whole point of the test.  Everything else
is the smallest configuration that still exercises ``read_topo_file``: one AMR
level, a few cells, one step.
"""

from clawpack.clawutil import data


def setrun(claw_pkg='geoclaw'):
    rundata = data.ClawRunData(claw_pkg, 2)

    cd = rundata.clawdata
    cd.num_dim = 2
    # Placeholders; the test overwrites lower/upper/num_cells before writing.
    cd.lower = [-95.0, 22.0]
    cd.upper = [-90.0, 25.0]
    cd.num_cells = [20, 12]
    cd.num_eqn = 3
    cd.num_aux = 3
    cd.capa_index = 2

    cd.t0 = 0.0
    cd.output_style = 1
    cd.num_output_times = 1
    cd.tfinal = 1.0
    cd.output_t0 = False
    cd.dt_initial = 0.1
    cd.cfl_desired = 0.75
    cd.steps_max = 100
    cd.verbosity = 1
    cd.dt_variable = True
    cd.bc_lower = ['extrap', 'extrap']
    cd.bc_upper = ['extrap', 'extrap']

    # The GeoClaw Riemann solver always returns 3 waves and uses f-waves.
    # These are not defaults: num_waves defaults to 1, so leaving them out lets
    # the solver write past the end of the wave arrays.  That corrupts memory
    # and shows up as an intermittent segfault *after* the run finishes -- with
    # fort.geo truncated, which looks like a topo problem and is not one.
    cd.order = 2
    cd.transverse_waves = 2
    cd.num_waves = 3
    cd.limiter = ['mc', 'mc', 'mc']
    cd.use_fwaves = True
    cd.source_split = 'godunov'

    amrdata = rundata.amrdata
    amrdata.amr_levels_max = 1
    amrdata.refinement_ratios_x = [2]
    amrdata.refinement_ratios_y = [2]
    amrdata.refinement_ratios_t = [2]
    amrdata.aux_type = ['center', 'capacity', 'yleft']

    geo = rundata.geo_data
    geo.gravity = 9.81
    geo.coordinate_system = 2
    geo.earth_radius = 6367.5e3
    geo.sea_level = 0.0
    geo.dry_tolerance = 1.e-3
    geo.friction_forcing = False

    rundata.refinement_data.wave_tolerance = 1.e-1

    return rundata


if __name__ == '__main__':
    setrun().write()
