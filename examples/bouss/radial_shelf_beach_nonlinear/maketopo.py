
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import Rearth  # radius of earth
from clawpack.clawutil.data import ClawData
from scipy.interpolate import interp1d

from numpy import *

from mapper import latlong, gcdist
x1d, z1d = loadtxt('1d_radial/celledges.data',skiprows=1,unpack=True)
z1d_func = interp1d(x1d, z1d, bounds_error=False, fill_value = z1d[-1])


def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 130*4 + 1
    nypoints = 130*4 + 1
    xlower= -130e3
    xupper=  130e3
    ylower= -130e3
    yupper=  130e3
    outfile= "ocean.topotype2"
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def makeqinit():
    """
    Create qinit data file
    """
    nxpoints=101
    nypoints=101
    xlower=-126.3
    xupper= -125.7
    ylower= 46.6
    yupper= 47.2
    outfile= "hump.xyz"
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def topo(x,y):
    """
    Cartesian: x,y in meters
    """
    import numpy as np
    x0 = 0.
    y0 = 0.
    #d = gcdist(x0,y0,x,y,Rearth)
    d = np.sqrt((x-x0)**2 + (y-y0)**2)
    z = z1d_func(d)
    return z



if __name__=='__main__':
    maketopo()
