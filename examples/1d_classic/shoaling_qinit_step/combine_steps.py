
from __future__ import print_function
from pylab import *
from scipy.interpolate import interp1d
import os

savefig_ext = '.png'
figdir = './figures'
os.system('mkdir -p %s' % figdir)

def save_figure(fname):
    """Save figure to figdir with desired extension"""
    full_fname = os.path.join(figdir,fname) + savefig_ext
    savefig(full_fname, bbox_inches='tight')
    print('Created %s' % full_fname)

outdir = '_output'
final_frame = 75
dp = loadtxt(outdir + '/fort.q%s' % str(final_frame).zfill(4), skiprows=5)
#dn = loadtxt('eta_neg.txt', skiprows=5)
g = loadtxt(outdir + '/celledges.data',skiprows=1)
xg = (g[:-1,0]+g[1:,0])/2.

eta_pos = dp[:,2]
#eta_neg = dn[:,2]
eta0 = interp1d(xg,eta_pos,fill_value="extrapolate")

def make_plot(frameno):
    figure(5, figsize=(7,6))
    clf()
    dn = loadtxt(outdir + '/fort.q%s' % str(frameno).zfill(4), skiprows=5)
    eta_neg = -dn[:,2]
    eta1 = interp1d(xg,eta_neg,fill_value="extrapolate")
    xx = linspace(-300e3,50e3,4000)
    subplot(211)
    plot(xx/1000,eta0(xx),'b')
    plot(xx/1000,eta1(xx),'g')
    title('Solution with positive step (blue) and negative step (green)')
    grid(True)
    xticks()
    ylabel('meters')
    subplot(212)
    plot(xx/1000, eta0(xx)+eta1(xx), 'r')
    title('Solution with square pulse (sum of step solutions)')
    grid(True)
    xlabel('kilometers')
    ylabel('meters')
    tight_layout()
    fname = 'steps_%s' % str(frameno).zfill(4)
    save_figure(fname)

if __name__=='__main__':
    make_plot(60)
    make_plot(70)
