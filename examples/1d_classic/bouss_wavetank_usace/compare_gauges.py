"""
Plot results at three gauges and compare with experimental data, obtained from
    https://nctr.pmel.noaa.gov/benchmark/Solitary_wave/

If you run this after running GeoClaw using the Madsen-Sorenson option with:
    rundata.bouss_data.boussEquations = 1
then you can compare the figure generated to Figure 4 in the BoussClaw paper
    http://dx.doi.org/10.1016/j.coastaleng.2017.01.005


See also 
    compare_BoussSWE.py
which runs the code with various settings and plots the results together
using plot_gauges from this module.

"""

from pylab import *
import clawpack.pyclaw.gauges as gauges

def plot_gauges(outdirs=[('_output', 'Simulation', 'b')],
                fname_figure='GaugeComparison.png'):

    d = loadtxt('experimental_data/ts3b.txt',skiprows=6)
    tg = d[:,0]
    g5 = d[:,2]
    g7 = d[:,4]
    g8 = d[:,5]

    figure(400, figsize=(8,8))
    clf()

    subplot(311)
    gaugeno = 5

    plot(tg-tg[0], g5, 'r--', label='Experiment')
    
    for (outdir,outdir_label,outdir_color) in outdirs:
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        t = gauge.t
        eta = gauge.q[2,:]
        plot(t, eta, color=outdir_color, label=outdir_label)


    xlim(0,25)
    ylim(-0.01,0.08)
    legend(loc='upper right')
    grid(True)
    xlabel('')
    ylabel('Surface (m)')
    title('Gauge %i' % gaugeno)

    subplot(312)
    gaugeno = 7
    
    plot(tg-tg[0], g7, 'r--', label='Experiment')
    
    for (outdir,outdir_label,outdir_color) in outdirs:
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        t = gauge.t
        eta = gauge.q[2,:]
        plot(t, eta, color=outdir_color, label=outdir_label)
        

    xlim(0,25)
    ylim(-0.01,0.08)
    legend(loc='upper right')
    grid(True)
    xlabel('')
    ylabel('Surface (m)')
    title('Gauge %i' % gaugeno)

    subplot(313)
    gaugeno = 8

    plot(tg-tg[0], g8, 'r--', label='Experiment')
    
    for (outdir,outdir_label,outdir_color) in outdirs:
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        t = gauge.t
        eta = gauge.q[2,:]
        plot(t, eta, color=outdir_color, label=outdir_label)


    xlim(0,25)
    ylim(-0.01,0.08)
    legend(loc='upper right')
    grid(True)
    xlabel('')
    ylabel('Surface (m)')
    title('Gauge %i' % gaugeno)

    tight_layout()

    if fname_figure is not None:
        savefig(fname_figure, bbox_inches='tight')
        print('Created %s' % fname_figure)

if __name__=='__main__':

    plot_gauges(outdirs=[('_output', 'Simulation', 'b')],
                fname_figure='GaugeComparison.png')

