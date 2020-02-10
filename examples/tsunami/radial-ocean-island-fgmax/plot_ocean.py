
from __future__ import print_function
from pylab import *

from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import Rearth  # radius of earth
from clawpack.clawutil.data import ClawData
from mapper import latlong
import maketopo

probdata = ClawData()
probdata.read('setprob.data', force=True)
theta_island = probdata.theta_island
print("theta_island = ",theta_island)


#close(1)
#figure(1, figsize=(10,6))
#axes([.1,.1,.5,.8])
clf()
s = linspace(0, 360., 1001)
xshore,yshore = latlong(1645.e3, s, 40., Rearth)
plot(xshore,yshore,'k',linewidth=2)

xshelf,yshelf = latlong(1560.e3, s, 40., Rearth)
plot(xshelf,yshelf,'k--',linewidth=2)

theta_island = 220.
(xi,yi) = latlong(1600.e3, theta_island, 40., Rearth)
xb = [xi-1, xi+1, xi+1, xi-1, xi-1]
yb = [yi-1, yi-1, yi+1, yi+1, yi-1]
plot(xb,yb,'b',linewidth=1.5)
text(5,49,'Test 1', fontsize=20)

theta_island = 260.
(xi,yi) = latlong(1600.e3, theta_island, 40., Rearth)
xb = [xi-1, xi+1, xi+1, xi-1, xi-1]
yb = [yi-1, yi-1, yi+1, yi+1, yi-1]
plot(xb,yb,'b',linewidth=1.5)
text(9,40,'Test 2', fontsize=20)

title('(a) Radial ocean',fontsize=20)

x = linspace(-5,5,51)
y = linspace(35,45,51)
X,Y = meshgrid(x,y)
Z = maketopo.qinit(X,Y)
contour(X,Y,Z,[2],linewidth=2,colors='k')

axis('scaled')
xlim([-21,21])
ylim([20,60])
xticks(fontsize='15')
yticks([25,35,45,55],fontsize='15')
xlabel('Longitude',fontsize=15)
ylabel('Latitude',fontsize=15)

savefig('ocean.png')
print("Created ocean.png")

if 0:
    axes([.65,.6,.25,.25])
    plot(xshore,yshore,'k',linewidth=2)
    plot(xshelf,yshelf,'k--',linewidth=2)
    maketopo.theta_island = 220.
    (xisland,yisland) = latlong(1600.e3, maketopo.theta_island, 40., Rearth)
    x = linspace(xisland-1, xisland+1, 200)
    y = linspace(yisland-1, yisland+1, 200)
    X,Y = meshgrid(x,y)
    Z = maketopo.topo(X,Y)
    contour(X,Y,Z,[0,7,15])
    axis('scaled')
    xlim([xisland-0.5, xisland+0.5])
    ylim([yisland-0.5, yisland+0.5])


    axes([.65,.2,.25,.25])
    plot(xshore,yshore,'k',linewidth=2)
    plot(xshelf,yshelf,'k--',linewidth=2)
    maketopo.theta_island = 260.
    (xisland,yisland) = latlong(1600.e3, maketopo.theta_island, 40., Rearth)
    x = linspace(xisland-1, xisland+1, 200)
    y = linspace(yisland-1, yisland+1, 200)
    X,Y = meshgrid(x,y)
    Z = maketopo.topo(X,Y)
    contour(X,Y,Z,[0,7,15])
    axis('scaled')
    xlim([xisland-0.5, xisland+0.5])
    ylim([yisland-0.5, yisland+0.5])

