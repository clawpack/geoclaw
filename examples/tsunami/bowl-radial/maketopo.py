
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 201
    nypoints = 201
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    outfile= "bowl.topotype2"     
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 101
    nypoints = 101
    xlower = -50.e0
    xupper = 50.e0
    yupper = 50.e0
    ylower = -50.e0
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    """
    Parabolic bowl
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    zmin = 80.
    z = 1.e-2*(x**2 + y**2) - zmin
    return z


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    ze = -((x+0e0)**2 + (y+0e0)**2)/10.
    z = where(ze>-10., 40.e0*exp(ze), 0.)
    return z

if __name__=='__main__':
    maketopo()
    makeqinit()
