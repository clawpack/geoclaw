"""
Plot fgmax output from GeoClaw run.

"""

import matplotlib.pyplot as plt
import os
import numpy
from clawpack.geoclaw import fgmax_tools
from clawpack.visclaw import geoplot

def plot_fgmax_grid(outdir,plotdir):

    fg = fgmax_tools.FGmaxGrid()
    fg.outdir = outdir
    data_file = os.path.join(outdir, 'fgmax_grids.data')
    fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
    fg.read_output()

    clines_zeta = [0.01] + list(numpy.linspace(0.05,0.3,6)) + [0.5,1.0,10.0]
    colors = geoplot.discrete_cmap_1(clines_zeta)
    plt.figure(1)
    plt.clf()
    zeta = numpy.where(fg.B>0, fg.h, fg.h+fg.B)   # surface elevation in ocean
    plt.contourf(fg.X,fg.Y,zeta,clines_zeta,colors=colors)
    plt.colorbar()
    plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k')  # coastline

    # plot arrival time contours and label:
    arrival_t = fg.arrival_time/3600.  # arrival time in hours
    clines_t = numpy.linspace(0,8,17)  # hours
    clines_t_label = clines_t[2::2]  # which ones to label 
    clines_t_colors = ([.5,.5,.5],)
    con_t = plt.contour(fg.X,fg.Y,arrival_t, clines_t,colors=clines_t_colors) 
    plt.clabel(con_t, clines_t_label)

    # fix axes:
    plt.ticklabel_format(style='plain',useOffset=False)
    plt.xticks(rotation=20)
    plt.gca().set_aspect(1./numpy.cos(fg.Y.mean()*numpy.pi/180.))
    plt.title("Maximum amplitude / arrival times")

    if not os.path.isdir(plotdir): 
        os.mkdir(plotdir)
    fname = os.path.join(plotdir, "amplitude_times.png")
    plt.savefig(fname)
    print("Created ",fname)

if __name__=="__main__":
    plot_fgmax_grid(outdir='_output', plotdir='.')
