"""
Plot timing info found in fort.t files
Requires modified valout function to print this info.
"""


from __future__ import print_function
from pylab import *
import glob

# Location of fort.t files:
outdir = '_output'

make_pngs = True  # print plots?

def make_png(fname):
    savefig(fname)
    #savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)

# set desired units for simulation time and computer time,
# based on length of run:

simtime_units = 'seconds'
comptime_units = 'seconds'

if simtime_units == 'seconds':
    simtime_factor = 1
elif simtime_units == 'minutes':
    simtime_factor = 60.
elif simtime_units == 'hours':
    simtime_factor = 3600.

if comptime_units == 'seconds':
    comptime_factor = 1
elif comptime_units == 'minutes':
    comptime_factor = 60.
elif comptime_units == 'hours':
    comptime_factor = 3600.


tfiles = glob.glob(outdir + '/fort.t*')
ntimes = len(tfiles)
print('Found %i fort.t files in %s' % (ntimes, outdir))

max_levels = 10  # never use more than this many levels

time = zeros(ntimes)
total_cpu = zeros(ntimes)
wtime = zeros((ntimes, max_levels))
cpu = zeros((ntimes, max_levels))
cells = zeros((ntimes, max_levels))

for j,fname in enumerate(tfiles):
    lines = open(fname).readlines()
    time[j] = float(lines[0].split()[0])
    try:
        total_cpu[j] = float(lines[11].split()[0])
    except:
        error_msg = '*** fort.t files do have have expected timing info'
        raise Exception(error_msg)
    #print(j, 'total cpu: ',total_cpu[j])
    for level in range(max_levels):
        lineno = 14+level
        if len(lines) < lineno+1:
            break
        #print(j, level, lines[lineno])
        tokens = lines[lineno].split()
        if len(tokens)==0:
            break
        wtime[j, level] = float(tokens[1])
        cpu[j, level] = float(tokens[2])
        cells[j, level] = float(tokens[3])
    
xlimits = [time.min()/simtime_factor, time.max()/simtime_factor]
ylimits = [0, 1.1*total_cpu.max()/comptime_factor]

if 0:
    figure(21)
    clf()
    plot(time/simtime_factor, total_cpu/comptime_factor)
    title('Total CPU time')
    xlabel('Simulation time (%s)' % simtime_units)
    ylabel('CPU time (%s)' % comptime_units)

figure(22)
clf()
sum_cells_over_levels = zeros(ntimes)
for j in range(max_levels):
    if max(cells[:,j]) == 0:
        break
    #plot(time/3600, cells[:,j], label='Level %s' % (j+1))
    last_sum_cells = sum_cells_over_levels.copy()
    sum_cells_over_levels += cells[:,j]
    fill_between(time/simtime_factor, last_sum_cells, sum_cells_over_levels, 
                 label='Level %s' % (j+1))

xlim(xlimits)
ylim(0, 1.1*sum_cells_over_levels[-1])
title('Cells updated on each level')
xlabel('Simulation time (%s)' % simtime_units)
ylabel('Grid cell updates')
legend(loc='upper left')

if make_pngs:
    make_png('CellUpdates.png')


figure(24)
clf()
sum_cpu_over_levels = zeros(ntimes)
for j in range(max_levels):
    if max(cpu[:,j]) == 0:
        break
    #plot(time/3600, cpu[:,j], label='Level %s' % (j+1))
    last_sum_cpu = sum_cpu_over_levels.copy()
    sum_cpu_over_levels += cpu[:,j]
    fill_between(time/simtime_factor, last_sum_cpu/comptime_factor, 
                 sum_cpu_over_levels/comptime_factor, 
                 label='Level %s' % (j+1))

plot(time/simtime_factor, total_cpu/comptime_factor, 'k', label='Total CPU')
xlim(xlimits)
ylim(ylimits)
title('CPU time on each level')
xlabel('Simulation time (%s)' % simtime_units)
ylabel('CPU time (%s)' % comptime_units)
legend(loc='upper left')

if make_pngs:
    make_png('CPUtime.png')


figure(25)
clf()
sum_wtime_over_levels = zeros(ntimes)
for j in range(max_levels):
    if max(wtime[:,j]) == 0:
        break
    last_sum_wtime = sum_wtime_over_levels.copy()
    sum_wtime_over_levels += wtime[:,j]
    fill_between(time/simtime_factor, last_sum_wtime/comptime_factor,
                 sum_wtime_over_levels/comptime_factor, 
                 label='Level %s' % (j+1))

title('Wall time on each level')
xlabel('Simulation time (%s)' % simtime_units)
ylabel('CPU time (%s)' % comptime_units)
legend(loc='upper left')
xlim(xlimits)
ylim(ylimits)

if make_pngs:
    make_png('WallTime.png')

