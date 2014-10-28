
import clawpack.geoclaw.shallow_1d.plot as geoplot

import setrun
rundata=setrun.setrun()


def setplot(plotdata):

    plotdata.clearfigures()

    plotfigure = plotdata.new_plotfigure(name='depth', figno=0)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0, 100.0]
    plotaxes.ylimits = [-1, 2]
    plotaxes.title = ''

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'

    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output

    return plotdata

