"""
Run GeoClaw with pure shallow water equations and with several choices
of Boussinesq equations, and then plot gauge results compared to experiment.

"""
from pylab import *
import os

from clawpack.clawutil.runclaw import runclaw
#from clawpack.visclaw.frametools import plotframe
#from clawpack.visclaw.data import ClawPlotData
import compare_gauges

import setrun

outdir_sgn = '_output_sgn'
print('outdir_bouss = ',outdir_sgn)

outdir_ms = '_output_ms'
print('outdir_ms = ',outdir_ms)

outdir_swe = '_output_swe'
print('outdir_swe = ',outdir_swe)

run_code = True  # set to False if output already exists

if run_code:
    # create executable and .data files:
    os.system('make .exe')
    rundata = setrun.setrun()

    # Boussinesq, MS:
    rundata.bouss_data.bouss_equations = 1
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_ms)   # run clawpack code

    # Boussinesq, SGN:
    rundata.bouss_data.bouss_equations = 2
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_sgn)   # run clawpack code
    
    # Shallow water equations:
    rundata.bouss_data.bouss_equations = 0
    rundata.write()
    runclaw(xclawcmd='xgeo',outdir=outdir_swe)   # run clawpack code
    

outdirs=[('_output_swe', 'SWE', 'k'), \
         ('_output_ms', 'MS', 'b'), \
         ('_output_sgn', 'SGN','g')]

compare_gauges.plot_gauges(outdirs, fname_figure='GaugeComparison_BoussSWE.png')
