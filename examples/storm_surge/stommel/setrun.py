# encoding: utf-8
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import clawpack.clawutil.oldclawdata as data
import numpy as np

# RAMP_UP_TIME = 0.0
RAMP_UP_TIME = 12 * 60**2

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

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    ndim = 2
    rundata = data.ClawRunData(claw_pkg, ndim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------

    rundata = setgeo(rundata)   # Defined below

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
    clawdata.ndim = ndim

    # Lower and upper edge of computational domain:
    clawdata.xlower = 0.0
    clawdata.xupper = 1000e3
    
    clawdata.ylower = 0.0
    clawdata.yupper = 1000e3

    # Number of grid cells:
    clawdata.mx = 100
    clawdata.my = 100

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 9

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0



    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.outstyle = 1
    
    # Number of days
    num_days = 100

    # Output interval per day
    step = 1.0
    
    if clawdata.outstyle==1:
        clawdata.nout = int(num_days / step) + int(np.ceil(clawdata.t0 / (step*60**2)))
        clawdata.tfinal = num_days * 24.0 * 60.0**2

    elif clawdata.outstyle == 2:
        # Specify a list of output times.
        # t_start = 0.7650e05
        # dt = 0.12E+03 
        # clawdata.tout = [0.0,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]
        # clawdata.tout = [x*60**2 for x in clawdata.tout]
        clawdata.tout = [0.0,1.0,2.0,3.0,4.0]
        # for i in xrange(60):
        #     clawdata.tout.append(i*dt+t_start)
        clawdata.nout = len(clawdata.tout)
    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 100
        clawdata.iout = [iout, ntot]



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 2



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 5000




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = 2

    # Number of waves in the Riemann solution:
    clawdata.mwaves = 3

    # List of limiters to use for each wave family:
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [3,3,3]

    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used,
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.src_split = 1


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.mthbc_xlower = 3
    clawdata.mthbc_xupper = 3

    clawdata.mthbc_ylower = 3
    clawdata.mthbc_yupper = 3


    # ---------------
    # AMR parameters:
    # ---------------


    # max number of refinement levels:
    mxnest = 1

    clawdata.mxnest = -mxnest   # negative ==> anisotropic refinement in x,y,t

    # List of refinement ratios at each level (length at least mxnest-1)
    # Run resolution.py 2 2 4 8 16 to see approximate resolutions
    clawdata.inratx = [2,2,2,2,2]
    clawdata.inraty = [2,2,2,2,2]
    clawdata.inratt = [2,2,2,2,2]
    # Instead of setting these ratios, set:
    # geodata.variable_dt_refinement_ratios = True
    # in setgeo.
    # to automatically choose refinement ratios in time based on estimate
    # of maximum wave speed on all grids at each level.


    # Specify type of each aux variable in clawdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    clawdata.auxtype = ['center','center','center','center','center','center',
                        'center','center','center']


    clawdata.tol = -1.0     # negative ==> don't use Richardson estimator
    clawdata.tolsp = 0.5    # used in default flag2refine subroutine
                            # (Not used in geoclaw!)

    clawdata.cutoff = 0.7   # efficiency cutoff for grid generation
    clawdata.kcheck = 3     # how often to regrid (every kcheck steps)
    clawdata.ibuff  = 2     # width of buffer zone around flagged points

    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geodata = rundata.geodata
    except:
        print "*** Error, this rundata has no geodata attribute"
        raise AttributeError("Missing geodata attribute")

    # == setgeo.data values ==
    geodata.variable_dt_refinement_ratios = True

    geodata.gravity = 9.81
    geodata.coordinate_system = 1
    geodata.earth_radius = 6367.5e3
    geodata.coriolis_forcing = True
    geodata.theta_0 = 45.0

    # == settsunami.data values ==
    geodata.dry_tolerance = 1.e-2
    geodata.wave_tolerance = 5.e-1
    geodata.speed_tolerance = [1e10,1e10,1e10,1e10,1e10]
    geodata.deep_depth = 2.e2
    geodata.max_level_deep = 4
    geodata.friction_forcing = 1
    geodata.manning_coefficient = 1e-6
    geodata.friction_depth = 1.e10

    # == settopo.data values ==
    geodata.topofiles = []
    # Analytic bathymetry
    geodata.topo_type = 1
    geodata.topo_location = 0.0
    geodata.topo_left = -1000.0
    geodata.topo_right = -1000.0
    # for topography, append lines of the form
    #   [topotype, minlevel, maxlevel, t1, t2, fname]

    # == setdtopo.data values ==
    geodata.dtopofiles = []
    # for moving topography, append lines of the form:  (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]
    # geodata.dtopofiles.append([1,3,3,'usgs100227.tt1'])

    # == setqinit.data values ==
    geodata.qinit_type = 0
    # geodata.qinit_type = 4
    geodata.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    # geodata.qinitfiles.append([1, 5, 'hump.xyz'])

    # == setregions.data values ==
    geodata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    # geodata.regions.append([1, 4, 0.0, 1e10, -99.00, -80.00, 17.00, 32.00]) # entire domain
    # geodata.regions.append([1, 5, 0.0, 1e10, -95.40, -94.42, 29.10, 29.92]) # Galveston Bay and inland
    # geodata.regions.append([1, 7, 0.0, 1e10, -94.84, -94.70, 29.30, 29.40]) # Channel into Galveston bay
    # geodata.regions.append([1, 7, 0.0, 1e10, -95.37, -95.9, 29.60, 29.83]) # Houston ship channel [-95º 22',-94º 54'] x [29º 36',29º 50']

    # == setgauges.data values ==
    geodata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]


    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]
    # geodata.fixedgrids.append([1e3,3.24e4,10,-90,-80,-30,-15,100,100,0,1])
    
    # == Multilayer ==
    geodata.layers = 1
    geodata.rho = 1000.0
    geodata.eta_init = 0.0
    geodata.richardson_tolerance = 0.95
    
    return rundata
    # end of function setgeo
    # ----------------------

def set_storm():

    import geoclaw.surge

    data = geoclaw.surge.StormData()

   # Physics parameters
    data.rho_air = 1.0              # Density of air (rho is not implemented above)
    data.ambient_pressure = 101.5e3 # Nominal atmos pressure

    # Source term controls
    data.wind_forcing = True
    data.pressure_forcing = True
    
    # Source term algorithm parameters
    data.wind_tolerance = 1e-6
    data.pressure_tolerance = 1e-4 # Pressure source term tolerance

    # AMR parameters, in m/s and meters respectively
    data.wind_refine = [20.0,40.0,60.0]
    data.R_refine = [60.0e3,40e3,20e3]
    
    # Storm parameters
    data.storm_type = 3 # Type of storm

    # Stommel wind field
    data.A = 0.2 / 1.2

    return data

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.write()

    # Create storm data
    storm_data = set_storm()
    storm_data.write()

