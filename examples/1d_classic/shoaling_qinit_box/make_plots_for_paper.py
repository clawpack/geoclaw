"""
Make figures 2 and 3 for journal article:

    Shoaling on Steep Continental Slopes: Relating Transmission and 
    Reflection Coefficients to Green's Law, 
    by J. George, D. I. Ketcheson, and R. J. LeVeque,
    Pure and Applied Geophysics, 2019.
    DOI 10.1007/s00024-019-02316-y

See http://faculty.washington.edu/rjl/pubs/Shoaling2019 for additional links.

"""
from pylab import *
import os

from clawpack.clawutil.runclaw import runclaw
from clawpack.visclaw.frametools import plotframe

from clawpack.visclaw.data import ClawPlotData
import make_celledges
from imp import reload
reload(make_celledges)

run_code = False  # set to False if output already exists
if run_code:
    # create executable and .data files:
    os.system('make .exe')
    os.system('make data')

savefig_ext = '.png'
figdir = './figures'
os.system('mkdir -p %s' % figdir)

def save_figure(fname):
    """Save figure to figdir with desired extension"""
    full_fname = os.path.join(figdir,fname) + savefig_ext
    savefig(full_fname, bbox_inches='tight')
    print('Created %s' % full_fname)

for xs in [0., 15e3, 60e3]:
    make_celledges.makegrid(xs)
    close('all')
    #import mapc2p
    import setplot
    #reload(mapc2p)
    reload(setplot)
    xs_str = 'box_xs%s' % str(int(xs/1e3)).zfill(2)
    w_str = 'w%s' % str(int(xs/1e3)).zfill(2)
    outdir = '_output_' + w_str
    print('outdir = ',outdir)
    if run_code:
        runclaw(xclawcmd='xgeo',outdir=outdir)   # run clawpack code
    pd = ClawPlotData()
    pd.outdir = outdir
    pd = setplot.setplot(pd)  # reload celledges.data for each xs value
    pd.outdir = os.path.abspath(outdir)
    pd.printfigs = False  # or else figure is closed after printing
    for frameno in [0,6]:
        plotframe(frameno,pd)
        fname = xs_str + '_frame%s' % str(frameno).zfill(2)
        save_figure(fname)


