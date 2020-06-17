"""
Plot fgmax output from GeoClaw run.

"""
try:
    matplotlib  # see if it's already been imported (interactive session)
except:
    import matplotlib
    matplotlib.use('Agg')  # set to image backend


from pylab import *
import matplotlib as mpl
import os
from clawpack.geoclaw import fgmax_tools


def add_gauges(label=True):
    fs = 15
    dy = -.002
    xy = [203.52825, 20.9021333]
    plot([xy[0]], [xy[1]], 'wo',markersize=8)
    plot([xy[0]], [xy[1]], 'k+',markersize=8)
    if label: text(203.529,20.9023,'1123',fontsize=15)

    xy = [203.530944, 20.895]
    plot([xy[0]], [xy[1]], 'wo',markersize=8)
    plot([xy[0]], [xy[1]], 'k+',markersize=8)
    if label: text(203.5293,20.8951,'TG',fontsize=15)

fg = fgmax_tools.FGmaxGrid()
# fg.read_input_data('fgmax1.txt')
fg.read_fgmax_grids_data(1)
fg.read_output(outdir='_output')

figure(1, figsize=(10,7))


bounds = 100*array([0,.25,.5,.75,1,2,4,5])  # cm/sec
cmap = mpl.colors.ListedColormap([[1,1,1],[.8,.8,1],[.5,.5,1],[0,0,1],\
                 [1,.7,.7], [1,.4,.4], [1,0,0]])
ax1 = axes()
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
contourf(fg.X,fg.Y,100*fg.s,bounds,cmap=cmap,norm=norm,extend='max')
cb = colorbar(extend='max')
cb.set_label('cm / sec')

contour(fg.X,fg.Y,100*fg.s,bounds,colors='w')
contour(fg.X,fg.Y,fg.B,[0],colors='k')

ticklabel_format(format='plain',useOffset=False)
#title('Maximum speed s')
xticks(rotation=20,fontsize=15)
yticks(fontsize=15)
ax1.set_aspect(1./cos(fg.Y.mean()*pi/180.))
xlim(203.515,203.5443)

add_gauges(False)

show()

# Plot gauges and topo:

bounds = [-1e10,0,1e10]
cmap = mpl.colors.ListedColormap([[1,1,1],[0,1,0]])
figure(2, figsize=(10,7))
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
contourf(fg.X,fg.Y,fg.B,bounds,cmap=cmap,norm=norm,extend='max')
contour(fg.X,fg.Y,fg.B,[0],colors='k')
ticklabel_format(format='plain',useOffset=False)
xticks(rotation=20,fontsize=15)
yticks(fontsize=15)
ax1.set_aspect(1./cos(fg.Y.mean()*pi/180.))
xlim(203.515,203.5443)

add_gauges()

show()

if 1:
    figure(2)
    fname = '../../Figures/Kahului_gauges.png'
    savefig(fname)
    print('Created ', fname)

    figure(1)
    fname = '../../Figures/Kahului_smax.png'
    savefig(fname)
    print('Created ', fname)
