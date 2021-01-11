
from __future__ import print_function
from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import Rearth  # radius of earth
from interp import pwcubic
from clawpack.clawutil.data import ClawData
import numpy as np
from matplotlib import pylab as plt

from mapper import latlong
import maketopo

probdata = ClawData()
probdata.read('setprob.data', force=True)
theta_island = probdata.theta_island
print("theta_island = ",theta_island)


(xisland,yisland) = latlong(1600.e3, theta_island, 40., Rearth)


def ocean():
    dmax = 1650.e3
    dp = np.linspace(0., dmax, 3301)
    zp = maketopo.shelf1(dp) 
    plt.clf()
    plt.plot(dp*1e-3,zp,'k-',linewidth=3)
    plt.fill_between(dp*1e-3,zp,0.,where=(zp<0.), color=[0,0.5,1])
    plt.xlim([0.,1650])
    plt.ylim([-4500.,500])
    plt.title("Topography as function of radius",fontsize=18)
    plt.xlabel("kilometers from center",fontsize=15)
    plt.ylabel("meters",fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig('topo1.png')
    #plt.savefig('topo1.tif')
    print("Cross-section plot in topo1.png")

    plt.xlim([1500,1650])
    plt.ylim([-200.,20])
    plt.title("Topography of shelf and beach",fontsize=18)
    plt.xlabel("kilometers from center",fontsize=15)
    plt.ylabel("meters",fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig('topo2.png')
    #plt.savefig('topo2.tif')
    print("Zoom near shore in topo2.png")


def with_island():
    dmax = 1650.e3
    dp = np.linspace(1500e3, dmax, 2000)
    zp = maketopo.shelf1(dp) + maketopo.island1(abs(dp-1600e3))

    plt.clf()
    plt.plot(dp*1e-3,zp,'k-',linewidth=3)
    plt.fill_between(dp*1e-3,zp,0.,where=(zp<0.), color=[0,0.5,1])

    plt.xlim([1500,1650])
    plt.ylim([-200.,30])
    plt.title("Cross-section through island center",fontsize=18)
    plt.xlabel("kilometers from center",fontsize=15)
    plt.ylabel("meters",fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig('topo3.png')
    #plt.savefig('topo3.tif')
    print("Cross-section through island in topo3.png")

if __name__=="__main__":
    ocean()
    with_island()
