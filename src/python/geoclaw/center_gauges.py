
from clawpack.geoclaw import center_points

def center_all_gaugedata(rundata, dx, dy, modify_rundata=True,
                         print_appends=True, verbose=False):
    """
    Function that can be called from setrun after setting rundata.gaugedata,
    to center all of the gauges using the domain edges and the specified dx,dy,
    the resolution at which you want them centered (in cells that will be
    an integer number of dx,dy away from the domain edges).
    
    Note that specifying any dx,dy less than the finest resolution in the
    simulation will at least insure that the gauges do not lie on grid
    cell edges on any AMR level (which can create noisy output if rounding
    errors cause values to be recorded at some times from one cell and at
    other times from its neighbor).
    """
    import copy

    if len(rundata.gaugedata.gauges) == 0:
        print('No gauges found to center')
        return rundata

    # assume the x,y values specified in gaugedata.gauges are the desired
    # locations that we will adjust a bit if necessary:
    # recall gaugedata.gauges is a list of lists for each gauge of the form:
    #    [gaugeno, x, y, t1, t2]

    x_desired = [gauge[1] for gauge in rundata.gaugedata.gauges]       
    y_desired = [gauge[2] for gauge in rundata.gaugedata.gauges]       
        
    # edges of cells should be integer multiples of dx,dy from domain edge:
    x_edge = rundata.clawdata.lower[0]
    y_edge = rundata.clawdata.lower[1]

    x_centered, y_centered = center_points.adjust_xy(x_desired, y_desired,
                                       x_edge, y_edge, dx, dy,
                                       verbose=verbose)

    # reset the x,y values for each gauge:
    gauges = copy.deepcopy(rundata.gaugedata.gauges)

    for k,gauge in enumerate(gauges):
        gauge[1] = float(x_centered[k])
        gauge[2] = float(y_centered[k])

    if print_appends:
        print('Modifications for setrun:')
        for k in range(len(gauges)):
            print('gauges.append([%i, %.8f, %.8f, %g, %g])' \
                    % tuple(gauges[k]))

    if modify_rundata:
        rundata.gaugedata.gauges = gauges
        print('*** Shifted rundata.gaugedata.gauges if necessary to center in cells')
    return rundata

if __name__=='__main__':

    # Sample code...

    # Center gauges based on rundata in setrun.py, assuming they 
    # should be centered at the finest AMR level specified by
    # `amr_level_max` and the refinement ratios.

    # Optional argument when executing for path to setrun file to use
    
    import sys
    from clawpack.clawutil.util import fullpath_import

    if len(sys.argv) > 1:
        setrun_file = sys.argv[1]
    else:
        setrun_file = 'setrun.py'

    print(f'Will center gauges at finest resolution based on setrun in {setrun_file}')
    setrun = fullpath_import(setrun_file)
        
        
    rundata = setrun.setrun()
    amr_levels_max = rundata.amrdata.amr_levels_max
    rrx = rundata.amrdata.refinement_ratios_x
    rry = rundata.amrdata.refinement_ratios_y
    lower = rundata.clawdata.lower
    upper = rundata.clawdata.upper
    dx_level1 = (upper[0] - lower[0]) / rundata.clawdata.num_cells[0]
    dy_level1 = (upper[1] - lower[1]) / rundata.clawdata.num_cells[1]

    dx_finest = dx_level1
    dy_finest = dy_level1
    for k in range(amr_levels_max-1):
        dx_finest = dx_finest / rrx[k]
        dy_finest = dy_finest / rry[k]

    print(f'Centering based on domain with:')
    print(f'          xlower = {lower[0]:.6f}, ylower = {lower[1]:.6f}')
    print(f'          dx_level1 = {dx_level1:.6e}, dy_level1 = {dy_level1:.6e}')
    print(f'          amr_levels_max = {amr_levels_max} with finest resolution')
    print(f'          dx_finest = {dx_finest:.6e}, dy_finest = {dy_finest:.6e}')

    # center gauges and print out the revised lines to insert in setrun.py
    center_all_gaugedata(rundata, dx_finest, dy_finest, modify_rundata=False,
                         print_appends=True, verbose=False)

    
