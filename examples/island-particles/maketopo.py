
"""
Module to create topo and qinit data files for this example.
"""

from pathlib import Path
import sys

import numpy as np

from clawpack.geoclaw.topotools import Topography

def maketopo(path=None):
    """
    Output topography file for the entire domain
    """

    if path:
        outfile = Path(path) / "island.tt3" 
    else:
        outfile = Path() / "island.tt3" 

    nxpoints = 201
    nypoints = 241
    xlower = 0.e0
    xupper = 100.e0
    ylower = 0.e0
    yupper = 50.e0

    topography = Topography(topo_func=topo)
    topography.x = np.linspace(xlower, xupper, nxpoints)
    topography.y = np.linspace(ylower, yupper, nypoints)
    topography.write(outfile, topo_type=3, Z_format="%22.15e")

def makeqinit(path=None):
    """
    Create qinit data file
    """
    if path:
        outfile = Path(path) / "qinit.xyz"
    else:
        outfile = Path() / "qinit.xyz"
    
    nxpoints = 101
    nypoints = 101
    xlower = -50.e0
    xupper = 50.e0
    yupper = 50.e0
    ylower = -50.e0

    topography = Topography(topo_func=qinit)
    topography.x = np.linspace(xlower, xupper, nxpoints)
    topography.y = np.linspace(ylower, yupper, nypoints)
    topography.write(outfile, topo_type=1)

def topo(x,y):
    """
    Island
    """
    ze = -((x-40.)**2 + (y-35.)**2)/20.
    
    #z_island = where(ze>-10., 100.*exp(ze), 0.)
    z_island = np.where(ze>-10., 150. * np.exp(ze), 0.)
    z = -50 + z_island
    return z


def qinit(x,y):
    """
    Dam break
    """
    return np.where(x<10, 40., 0.)

if __name__=='__main__':
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = Path()
    maketopo(path)
    makeqinit(path)
