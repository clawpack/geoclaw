"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

Leave this version alone for reading into the Jupyter notebook to illustrate
various settings.

"""

from __future__ import absolute_import
from __future__ import print_function
import os, sys
import numpy as np


try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

from clawpack.geoclaw.data import ForceDry
from clawpack.amrclaw.data import FlagRegion


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
    #probdata.add_param('variable_eta_init', True)  # now in qinit info


    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
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
    # x values should be integer multipes of 1/3"
    # y values should be integer multipes of 1/3"
    # Note: always satisfied if limits are multiples of 0.01 degree

    arcsec16 = 1./(6*3600.)

    # choose domain and offset edges by half a 1/3" cell so
    # cell centers are exactly at DEM grid points:

    clawdata.lower[0] = -1.9 - arcsec16     # west longitude
    clawdata.upper[0] = 0.1 - arcsec16     # east longitude

    clawdata.lower[1] = -1.9 - arcsec16      # south latitude
    clawdata.upper[1] = 1.9 - arcsec16        # north latitude

    # choose mx and my so coarsest grid has 2 minute resolution:
    clawdata.num_cells[0] = 60
    clawdata.num_cells[1] = 114


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
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False   # True to restart from prior results
    clawdata.restart_file = ''

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 5
        clawdata.tfinal = 25*60.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 20
        clawdata.output_t0 = True
        

    clawdata.output_format = 'binary'

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



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

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.8 
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
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

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
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'



    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # negative checkpoint_style means alternate between aaaaa and bbbbb files
    # so that at most 2 checkpoint files exist at any time, useful when
    # doing frequent checkpoints of large problems.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = 3600.*np.arange(1,16,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least mxnest-1)
    # dx = dy = 2', 10", 2", 1/3":
    amrdata.refinement_ratios_x = [12,5,6]
    amrdata.refinement_ratios_y = [12,5,6]
    amrdata.refinement_ratios_t = [12,5,6]



    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 1  


    # ---------------
    # Regions:
    # ---------------
    #rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    # NO OLD STYLE REGIONS USED HERE


    # ---------------
    # NEW flagregions
    # ---------------


    flagregions = rundata.flagregiondata.flagregions  # initialized to []

    # now append as many flagregions as desired to this list:

    # The entire domain restricted to level 1 for illustration:
    # Note that this is a rectangle specified in the new way:
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 1
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    # domain plus a bit so kml files look nicer:
    flagregion.spatial_region = [clawdata.lower[0] - 0.1,
                                 clawdata.upper[0] + 0.1,
                                 clawdata.lower[1] - 0.1,
                                 clawdata.upper[1] + 0.1]
    flagregions.append(flagregion)

    # force 2 levels around dtopo source region for short time:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_level2_dtopo'
    flagregion.minlevel = 2
    flagregion.maxlevel = 2
    flagregion.t1 = 0.
    flagregion.t2 = 2.
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-2,1,-1,1]
    flagregions.append(flagregion)

    # allow 3 levels around coastal region for all times:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_level3'
    flagregion.minlevel = 1
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-0.01,0.01,-0.01,0.01]
    flagregions.append(flagregion)

    # force 4 levels around coastal region starting at 5 minutes:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_level4'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 5*60.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-0.005, 0.01, -0.011, 0.011]
    flagregions.append(flagregion)

    # ---------------
    # Gauges:
    # ---------------
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    rundata.gaugedata.gauges = []


    # Set GeoClaw specific runtime parameters.

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.2

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]
    topodir = 'input_files'

    topofiles.append([3, topodir + '/topo_ocean.tt3'])
    topofiles.append([3, topodir + '/topo_shore.tt3'])


    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]
    dtopodir = 'input_files'
    dtopo_data.dtopofiles.append([3, dtopodir + '/dtopo_test.tt3'])

    dtopo_data.dt_max_dtopo = 1.0


    # == setqinit.data values ==
    
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [fname]
    
    # initially set to False, will be changed in RunGeoClaw.ipynb:
    rundata.qinit_data.variable_eta_init = False 
    
    if 0:
        # initially no force_dry array, added in RunGeoClaw.ipynb
        force_dry = ForceDry()
        force_dry.tend = 7*60.
        force_dry.fname = 'input_files/force_dry_init.tt3'
        rundata.qinit_data.force_dry_list.append(force_dry)

    # == fgmax.data values ==
    #fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files


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
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    return rundata
    # end of function setrun
    # ----------------------



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    

