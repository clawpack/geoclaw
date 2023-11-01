"""
Plot results at gauges to compare with Figure 5 of Matsuyama et al. (2007).
Blue curves are the GeoClaw results.

Wavetank data for some of the gauges is available in the file
data_wavegauge.csv (kindly provided by Prof. Matsuyama).
At these gauges the observations are also plotted as red curves.
"""

from pylab import *
import clawpack.pyclaw.gauges as gauges

outdir = '_output'

add_data = True


figure(400, figsize=(14,8))
clf()

xp_gauges = [-150, -80., -60, -50, -40, -32, -31.2, -30.8, -30, \
             -28, -20, -10, -5]
             
gauge_info = [(-80,   1,(10,60)),
              (-60,   3,(20,70)),
              (-50,   5,(30,80)),
              (-40,   7,(40,90)),
              (-32,   9,(40,90)),
              (-31.2,11,(40,90)),
              (-30.8, 2,(40,90)),
              (-30,   4,(40,90)),
              (-28,   6,(40,90)),
              (-20,   8,(50,100)),
              (-10,  10,(60,110)),
              (-5,   12,(60,110))]
          
                 
for info in gauge_info:
    xp_g, k, xlimits = info
    subplot(6,2,k)
    gaugeno = int(-xp_g*10)
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    t = gauge.t
    eta = gauge.q[2,:]

    plot(t, 100*eta, 'b', label='x = %.1f' % xp_g)

    ylim(-5, 15)
    xlim(xlimits)
    grid(True)
    yticks([-5,0,5,10,15])
    xlabel('')
    #ylabel('Surface (m)')
    #title('Gauge %i' % gaugeno)

    if k != 2:
        legend(loc='upper right',framealpha=1)

if add_data:
    d = loadtxt('data_wavegauge.csv',skiprows=2,delimiter=',')
    subplot(6,2,5)
    plot(d[:,0],d[:,2],'r')
    subplot(6,2,7)
    plot(d[:,0],d[:,3],'r')
    subplot(6,2,2)
    plot(d[:,0],d[:,4],'r',label='observed')
    legend(loc='upper right',framealpha=1)
    subplot(6,2,8)
    plot(d[:,0],d[:,5],'r')
    subplot(6,2,12)
    plot(d[:,0],d[:,6],'r')

tight_layout()

if 1:
    fname = 'GaugesFigure5.png'
    savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)
