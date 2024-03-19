
# OLD:

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('minLevelBouss', 1,' minlevel for Bouss terms')
    probdata.add_param('maxLevelBouss', 10,' maxlevel for Bouss terms')
    probdata.add_param('deepBouss', 5,' min water depth for Bouss terms')
    probdata.add_param('solver', 3,     ' 1=GMRES, 2=Pardiso, 3=PETSc')
    probdata.add_param('equations', 2,     ' 1=MadsenSchaffer, 2=SGN')
    probdata.add_param('alpha', 1.153,  ' If using SGN, else ignore')
    probdata.add_param('startWithBouss', True, 'Take numSWEsteps of SWE first')
    probdata.add_param('numSWEsteps', 0,  ' Take this many SWE steps first')


# NEW:

    # To use Boussinesq solver, add bouss_data parameters here
    # Also make sure to use the correct Makefile pointing to bouss version
    from clawpack.geoclaw.data import BoussData
    rundata.add_data(BoussData(),'bouss_data')
    
    # CHECK ORDER!

    rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN
    rundata.bouss_data.bouss_min_level = 1    # coarsest level to apply bouss
    rundata.bouss_data.bouss_max_level = 10   # finest level to apply bouss
    rundata.bouss_data.bouss_min_depth = 10.  # depth to switch to SWE
    rundata.bouss_data.bouss_solver = 3       # 1=GMRES, 2=Pardiso, 3=PETSc
    rundata.bouss_data.bouss_tstart = 0.      # time to switch from SWE
