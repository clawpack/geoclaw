"""
Make plots for Figures 5 and 6 in the paper.
"""

from __future__ import print_function
from pylab import *
from scipy.integrate import quad
import os
from clawpack.clawutil.runclaw import runclaw

import combine_steps  

run_code = True  # set to fault if output already exists
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


xs = 10.e3
xlimits = [-150e3,50e3]
outdir = '_output'

if run_code:
    runclaw(xclawcmd='xgeo',outdir=outdir)   # run clawpack code


griddata = loadtxt(outdir + '/celledges.data', skiprows=1)
xgrid = griddata[:,0]
zgrid = griddata[:,1]
xcell = 0.5*(xgrid[:-1]+xgrid[1:])

figure(11, figsize=(7,3))
clf()
#plot(xgrid,zgrid,'k')
if 1:
    fill_between(xgrid,zgrid,0,color=[0.7,0.7,1])
    plot(xgrid,zgrid,'k')
    plot(xgrid,0*zgrid,'b')
    ylim(-3500,300)
#title('Bathymetry')
xticks([-xs,0,xs],['$-\epsilon$','0','$\epsilon$'])
xlim(xlimits)
ylabel('meters')

save_figure('wave_topo')


framenos = [0,5,10,20,40,70]
for frameno in framenos:
    fname = outdir + '/fort.q%s' % str(frameno).zfill(4)
    q = loadtxt(fname, skiprows=6)
    fname = outdir + '/fort.t%s' % str(frameno).zfill(4)
    t = float(open(fname).readline().split()[0])
    print('t = %g' % t)
    figure(12, figsize=(7,3))
    clf()
    plot(xcell, q[:,2], 'b')

    txs = t/xs
    if t==0:
        title('t = 0')
    else:
        title('t = %4.3f$ \epsilon$' % txs)
    xticks([-xs,0,xs],['$-\epsilon$','0','$\epsilon$'])

    hl = 3200.
    hr = 200.
    greens = (hl/hr)**(0.25)
    #print 'greens = ',greens
    #plot(current_data.x, greens*ones(current_data.x.shape),'g--')
    plot(xlimits,[greens,greens],'g--',label="$C_G$, Green's Law")
    ctrans = 2*sqrt(hl)/(sqrt(hl)+sqrt(hr))
    plot(xlimits,[ctrans,ctrans],'r--',label="$C_T$, Transmission coefficient")
    #crefl = (sqrt(hl)-sqrt(hr))/(sqrt(hl)+sqrt(hr))
    #print 'ctrans = ',ctrans
    legend(loc='lower left')

    xlim(xlimits)
    ylabel('meters')
    draw()

    fname = 'wave_%s' % str(frameno).zfill(4)
    save_figure(fname)

# make plots for Figure 6 in paper:

combine_steps.make_plot(60)
combine_steps.make_plot(70)
