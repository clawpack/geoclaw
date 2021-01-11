""" 
Module to set up run time parameters for Clawpack -- AMRClaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

from __future__ import print_function
import os, sys
import numpy as np

from clawpack.geoclaw.data import Rearth  # radius of earth
from mapper import latlong

from clawpack.geoclaw import fgmax_tools

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

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    theta_island = 220.
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('theta_island', theta_island,  'angle to island')

    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -20.0          # xlower
    clawdata.upper[0] = 20.0          # xupper
    clawdata.lower[1] = 20.0          # ylower
    clawdata.upper[1] = 60.0          # yupper
    
    # Number of grid cells:
    clawdata.num_cells[0] = 40      # mx
    clawdata.num_cells[1] = 40      # my
    

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
 
    clawdata.output_style = 1
 
    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 7
        clawdata.tfinal = 14000.
        clawdata.output_t0 = True  # output at initial (or restart) time?
        
    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  [0., 0.1]
 
    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 2
        clawdata.total_steps = 4
        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    clawdata.output_format = 'binary'      # 'ascii', 'binary'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 16.0
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer']
    
    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper
                  
       
    # ---------------
    # Gauges:
    # ---------------
    gauges  = rundata.gaugedata.gauges 
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    gaugeno = 0
    for d in [1570e3, 1590e3, 1610e3, 1630e3]:
        gaugeno = gaugeno+1
        x,y = latlong(d, theta_island, 40., Rearth)
        gauges.append([gaugeno, x, y, 0., 1e10])

                  
    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 1

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.  
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5

    

    # ---------------
    # AMR parameters:   (written to amr.data)
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 5

    # List of refinement ratios at each level (length at least amr_level_max-1)
    amrdata.refinement_ratios_x = [4, 3, 5, 4]
    amrdata.refinement_ratios_y = [4, 3, 5, 4]
    amrdata.refinement_ratios_t = [1, 1, 1, 1]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'capacity', 'yleft']


    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance
    
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.5  # tolerance used in this routine
    # Note: in geoclaw the refinement tolerance is set as wave_tolerance below 
    # and flag2refine_tol is unused!

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3       

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0      


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    # default on full domain (no longer implicit in topo file, as of v5.8.0):
    regions.append([1, 1,    0., 1e9, -20., 20., 20., 60.])
    # initially around hump (no longer implicit in qinit file, as of v5.8.0):
    regions.append([1, 2,    0., 1., -20., 20., 20., 60.])

    regions.append([1, 3,    0., 5000., -5., 20., 35., 55.])
    regions.append([1, 2, 5000., 6900., -5., 20., 35., 55.])
    regions.append([1, 3, 5000., 6900., 10., 20., 45., 55.])
    regions.append([1, 2, 6900., 9000., 10., 20., 47., 52.])

    # Force refinement near the island as the wave approaches:

    (xisland,yisland) = latlong(1600.e3, theta_island, 40., Rearth)
    x1 = xisland - 1.
    x2 = xisland + 1.
    y1 = yisland - 1.
    y2 = yisland + 1.
    regions.append([4, 4, 7000., 1.e10,  x1,x2,y1,y2])

    x1 = xisland - 0.2
    x2 = xisland + 0.2
    y1 = yisland - 0.2
    y2 = yisland + 0.2
    regions.append([4, 5, 8000., 1.e10,  x1,x2,y1,y2])

    # -----------------------------------------------
    # GeoClaw specific parameters:

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system =  2
    geo_data.earth_radius = 6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 0.001
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1000000.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.01

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]
    topofiles.append([3, 'ocean.tt3'])
    topofiles.append([3, 'island.tt3'])

    # == setdtopo.data values ==
    rundata.dtopo_data.dtopofiles = []
    dtopofiles = rundata.dtopo_data.dtopofiles
    # for moving topography, append lines of the form :  
    #   [topotype, fname]


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  4
    qinitfiles = rundata.qinit_data.qinitfiles 
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [fname]
    rundata.qinit_data.qinitfiles = [['hump.xyz']]


    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2  # Save depth and speed

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.
    
    # Here several different ones are specified to illustrate:
    
    tstart_max = 8000.  # when to start monitoring fgmax grids

    # Points on a uniform 2d grid:
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid  
    fg.x1 = 14.25
    fg.x2 = 14.65
    fg.y1 = 50.10
    fg.y2 = 50.35
    fg.dx = 15/ 3600.  # desired resolution of fgmax grid
    fg.min_level_check = amrdata.amr_levels_max # which levels to monitor max on
    fg.tstart_max = tstart_max  # just before wave arrives
    fg.tend_max = 1.e10    # when to stop monitoring max values
    fg.dt_check = 20.      # how often to update max values
    fg.interp_method = 0   # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)  # written to fgmax_grids.data

    # Points on a 1d transect from (x1,y1) to (x2,y2):
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 1  # equally spaced points on a transect
    fg.x1 = 14.25
    fg.x2 = 14.65
    fg.y1 = 50.10
    fg.y2 = 50.35
    fg.npts = 50
    fg.min_level_check = amrdata.amr_levels_max # which levels to monitor max on
    fg.tstart_max = tstart_max  # just before wave arrives
    fg.tend_max = 1.e10    # when to stop monitoring max values
    fg.dt_check = 20.      # how often to update max values
    fg.interp_method = 0   # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)  # written to fgmax_grids.data

    # fgmax grid point_style==4 means grid specified as topo_type==3 file:
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 4
    fg.min_level_check = amrdata.amr_levels_max # which levels to monitor max on
    fg.tstart_max = tstart_max  # just before wave arrives
    fg.tend_max = 1.e10    # when to stop monitoring max values
    fg.dt_check = 20.      # how often to update max values
    fg.interp_method = 0   # 0 ==> pw const in cells, recommended
    fg.xy_fname = 'fgmax_pts_island.data'  # file of 0/1 values in tt3 format
    fgmax_grids.append(fg)  # written to fgmax_grids.data

    # fgmax grid point_style==0 means list of points:
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 0
    fg.min_level_check = amrdata.amr_levels_max # which levels to monitor max on
    fg.tstart_max = tstart_max  # just before wave arrives
    fg.tend_max = 1.e10    # when to stop monitoring max values
    fg.dt_check = 20.      # how often to update max values
    fg.interp_method = 0   # 0 ==> pw const in cells, recommended
    # can set list of points here:
    fg.npts = 2
    fg.X = np.array([14.4, 14.5])
    fg.Y = np.array([50.13, 50.13])
    fgmax_grids.append(fg)  # written to fgmax_grids.data

    # fgmax grid point_style==0 means list of points:
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 0
    fg.min_level_check = amrdata.amr_levels_max # which levels to monitor max on
    fg.tstart_max = tstart_max  # just before wave arrives
    fg.tend_max = 1.e10    # when to stop monitoring max values
    fg.dt_check = 20.      # how often to update max values
    fg.interp_method = 0   # 0 ==> pw const in cells, recommended
    # can specify that list of points is in a different file:
    fg.npts = 0
    fg.xy_fname = 'fgmax_ps0.txt'
    fgmax_grids.append(fg)  # written to fgmax_grids.data
    
    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    
    # remake topo in case island location theta_island was changed above:
    #import maketopo
    #maketopo.maketopo()
    #maketopo.makeqinit()
    
