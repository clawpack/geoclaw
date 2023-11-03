"""
Make plots for the paper.
Note that figno=120 is used, as specified in setplot.py
and the xlimits,ylimits for the specific frames used here are set there.
The savefig to a pdf file is also done in setplot.
"""

from pylab import *
import setplot

plotdata = setplot.setplot()

plotdata.printfigs = False
plotdata.print_fignos = [120]

plotdata.outdir = '_output'

for frameno in [6,12,18]:
    plotdata.plotframe(frameno)

# For later times,
# reset plotitems in pcolor plot so coarser grids aren't shown near shore:
plotfigure = plotdata.plotfigure_dict['For paper']
plotaxes = plotfigure.plotaxes_dict['pcolor']
for plotitem in plotaxes.plotitem_dict.values():
    plotitem.amr_data_show = [0,0,1,1]
    
for frameno in [21,23,26]:
    plotdata.plotframe(frameno)
