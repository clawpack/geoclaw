"""
Module to set up run time parameters for geoclaw 1d_nonuniform code

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os, sys
import numpy as np
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p


# Read in nonuniform computational cell edges, which should have
# been created using make_celledges.py:

grid_type = 2
fname_celledges = 'celledges.data'
    
mapc2p, mx_edge, xp_edge = make_mapc2p(fname_celledges)
mx = mx_edge - 1
        
print('Setting mx = %i, cell edges from %g to %g' % (mx,xp_edge[0],xp_edge[-1]))


#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data


    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 1
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    # Sample setup to write one line to setprob.data ...
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')


    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    # For nonuniform grid, 0 <= xc <= 1 and the file grid.data should
    # define the mapping to the physical domain

    clawdata.lower[0] = 0.          # xlower
    clawdata.upper[0] = 1.           # xupper

    # Number of grid cells:
    clawdata.num_cells[0] = mx
    
    # 1D geoclaw requires grid_data:
    rundata.grid_data.grid_type = grid_type  # should be set to 2 above
    rundata.grid_data.fname_celledges = fname_celledges
    rundata.grid_data.monitor_fgmax = False  # record max h,s,etc in each cell?
    rundata.grid_data.monitor_runup = False  # record first and last wet cells?
    rundata.grid_data.monitor_total_zeta = False # record "total mass in wave"?

    # To use Boussinesq solver, add bouss_data parameters here
    # Also make sure to use the correct Makefile pointing to bouss version
    from clawpack.geoclaw.data import BoussData1D
    rundata.add_data(BoussData1D(),'bouss_data')

    rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN
    rundata.bouss_data.bouss_min_depth = 0.06  # depth to switch to SWE

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 2

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 2

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2


    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.qNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.q0006'   # File to use for restart data


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 48
        clawdata.tfinal = 120.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  []

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 10
        clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format = 'ascii'  # 'ascii' is only option currently

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'all'  # could be list
    clawdata.output_aux_onlyonce = True   # output aux arrays only at t0


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==False: fixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.01

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.e9

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.25
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2


    # Number of waves in the Riemann solution:
    clawdata.num_waves = 2

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = [4,4]

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bc1.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'user'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper


    geo_data = rundata.geo_data

    geo_data.dry_tolerance = 1.e-6

    # Friction source terms:
    #   src_split > 0 required
    #   currently only Manning friction with a single n=friction_coefficient
    #   is supported in 1d.

    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025

    geo_data.coordinate_system = 1  # linear distance (meters)

    topo_data = rundata.topo_data
    topo_data.topofiles.append([1, 'celledges.data'])


    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, t1, t2]

    # for gauges append [gauge id, xc, t1, t2])
    # note that xc is the computational grid point, 0 <= xc <= 1,
    # so if you want to specify physical points xp, these need to be mapped
    # to corresponding xc as follows:

    if 1:
        #xp_gauges = [-150, -80., -40, -30, -20, -10, 0]
        xp_gauges = [-150, -80., -60, -50, -40, -32, -31.2, -30.8, -30, \
                     -28, -20, -10, -5]
        for k,xp_g in enumerate(xp_gauges):
            #gaugeno = k+1  
            gaugeno = int(-xp_g*10)
            # compute computational point xc_g that maps to xp_g:
            ii = np.where(xp_edge < xp_g)[0][-1]
            xp_frac = (xp_g - xp_edge[ii])/(xp_edge[ii+1] - xp_edge[ii])
            xc_g = (ii + xp_frac)/float(mx)
            print('gaugeno = %i: physical location xp_g = %g maps to xc_g = %.12f' \
                  % (gaugeno,xp_g, xc_g))
            rundata.gaugedata.gauges.append([gaugeno, xc_g, 0, 1e9])



    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

