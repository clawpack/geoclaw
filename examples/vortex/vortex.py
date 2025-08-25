#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: flow over a sill
==================================

Solve the 2D shallow water equations with
variable bathymetry:

.. :math:
    h_t + (hu)_x + (hv)_y & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y & = -g h b_x \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y & = -g h b_y.

The bathymetry is flat in this example, but there is a translating vortex
"""

import sys
import os
import subprocess

import numpy as np

from clawpack import riemann
from clawpack import pyclaw
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn

import setplot

def setup(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          outdir='./_py_output', N=256):

    solver = pyclaw.ClawSolver2D(riemann.shallow_bathymetry_fwave_2D)
    solver.dimensional_split = 1  # No transverse solver available

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.extrap
    solver.aux_bc_upper[1] = pyclaw.BC.extrap

    x = pyclaw.Dimension(-50, 50, N, name='x')
    y = pyclaw.Dimension(-50, 50, N, name='y')
    domain = pyclaw.Domain([x, y])
    state = pyclaw.State(domain, num_eqn, num_aux=1)

    X, Y = state.p_centers
    state.aux[0,:,:] = -np.ones((N, N))

    h, u, v = setplot.exact_solution(X, Y, 0.0)
    state.q[depth,:,:] = h
    state.q[x_momentum,:,:] = h * u
    state.q[y_momentum,:,:] = h * v

    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1.e-3
    state.problem_data['sea_level'] = 0.

    claw = pyclaw.Controller()
    claw.tfinal = 100
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.setplot = setplot
    claw.outdir = outdir
    claw.write_aux_init = True
    claw.keep_copy = True

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main, _info_from_argv
    petsc_args, pyclaw_kwargs = _info_from_argv(sys.argv)
    run_conv = bool(pyclaw_kwargs.pop('conv', False))
    if run_conv:
        base_path = os.path.expandvars(os.path.join("$DATA_PATH", "vortex_example", "pyclaw"))
        num_cells = [2**n for n in range(6, 11)]
        for N in num_cells:
            path = os.path.join(base_path, f"n{str(N).zfill(4)}_output")
            claw = setup(outdir=path, N=N)
            print(f"Running {N}: {path}")
            claw.run()
        subprocess.run(["./plot_comparison.py", num_cells])
    else:
        claw = run_app_from_main(setup, setplot)

    
