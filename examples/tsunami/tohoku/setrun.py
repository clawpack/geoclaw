"""
Module to set up run time parameters for Clawpack -- AMRClaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW environment variable")

topodir = os.path.join(CLAW, 'geoclaw', 'scratch')
if not os.path.isdir(topodir):
    raise Exception("*** Missing topo directory: %s" % topodir)


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
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')


    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

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
    clawdata.lower[0] = 132.0          # xlower
    clawdata.upper[0] = 210.0          # xupper
    clawdata.lower[1] = 9.0          # ylower
    clawdata.upper[1] = 53.0          # yupper

    # Number of grid cells:
    clawdata.num_cells[0] = 39      # mx
    clawdata.num_cells[1] = 22      # my


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
        clawdata.num_output_times = 26
        clawdata.tfinal = 13*3600.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  np.linspace(7,13,4)*3600.

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = False  # output at initial (or restart) time?


    clawdata.output_format == 'binary'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0


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
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000


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
    gauges = rundata.gaugedata.gauges

    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    gauges.append([1123, 203.52825, 20.9021333, 7.0*3600., 1.e9]) #Kahului

    # more accurate coordinates from Yong Wei at PMEL:
    gauges.append([5680, 203.530944, 20.895, 7.0*3600., 1.e9]) #TG Kahului


    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

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
    amrdata.amr_levels_max = 5   # Set to 6 originally

    # List of refinement ratios at each level (length at least amr_level_max-1)
    # 2 degree, 24', 4', 1', 10", 1/3"
    amrdata.refinement_ratios_x = [5, 6, 4, 6, 30]
    amrdata.refinement_ratios_y = [5, 6, 4, 6, 30]
    amrdata.refinement_ratios_t = [5, 6, 4, 6, 30]


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

    inf = 1e9

    # Region 0 : Global region.  This assures a maximum refinement in any 
    # regions not covered by other regions listed below.
    regions.append([1, 2, 0., inf, 0, 360, -90, 90])

    # Region 1 : (from the dtopo file, below).  
    # Time interval taken from dtopo file : [0,1]
    regions.append([3, 3, 0, 1, 135., 150., 30., 45.])

    # Region 2 : Large region encompassing lower 3/4 of domain.
    # Time interval :  (0,18000)
    regions.append([1, 3, 0., 5.*3600., 132., 220., 5., 40.])

    # Region 3 : Region including all Hawaiian Islands.
    # Time interval  (18000,28800)
    regions.append([1, 3, 5.*3600.,  8.*3600., 180., 220., 5., 40.])

    # Region 4 : Region including Molekai and Maui. 
    # Time interval : (23400.0, inf)
    regions.append([4, 4, 6.5*3600., inf, 202.5,204,20.4,21.4])

    # Region 5 : Strip including north shore of Maui.
    # Time interval :  (25200, inf)
    regions.append([5, 5, 7.*3600., inf, 203.0, 203.7, 20.88333, 21.])

    # Region 6 :  Includes port at Kailua.
    # Time interval :  (26100.0, inf)
    regions.append([6, 6,  7.25*3600., inf, 203.52,    203.537, 20.89,   20.905])


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


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    """

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
    geo_data.manning_coefficient = 0.035
    geo_data.friction_depth = 500.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.02
    refinement_data.deep_depth = 200.0
    refinement_data.max_level_deep = 4

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    topofiles.append([3, 1, 1, 0.0, 1e10, os.path.join(topodir,'etopo1min130E210E0N60N.asc')])
    # topofiles.append([3, 1, 1, 0.0, 1e10, os.path.join(topodir,'hawaii_6s.txt')])
    topofiles.append([3, 1, 1, 0., 1.e10, os.path.join(topodir,'kahului_1s.txt')])


    # == setdtopo.data values ==
    # [dtype, minlevel, maxlevel, 'topo']

    # Region 1, above handles this region.  
    rundata.dtopo_data.dtopofiles = [[1, 1, 1, os.path.join(topodir,'fujii.txydz')]]

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  0
    rundata.qinit_data.qinitfiles = []

    # == fixedgrids.data values ==
    # fixedgrids = []
    # rundata.fixed_grid_data.fixedgrids = []
    # fixedgrids = rundata.fixed_grid_data.fixedgrids

    # == fgmax.data values ==
    # fgmax_files = rundata.fgmax_data.fgmax_grids
    # for fixed grids append to this list names of any fgmax input files
    # fgmax_files.append('fgmax1.txt')
    # rundata.fgmax_data.num_fgmax_val = 2


    return rundata
    # end of function setgeo
    # ----------------------

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools
    rundata = setrun(*sys.argv[1:])
    rundata.write()

    kmltools.make_input_data_kmls(rundata)
