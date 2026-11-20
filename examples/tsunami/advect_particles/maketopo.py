
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
        outfile = Path(path) / "flat.tt3" 
    else:
        outfile = Path() / "flat.tt3" 

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
    

    nxpoints = 201
    nypoints = 241
    xlower = 0.e0
    xupper = 100.e0
    ylower = 0.e0
    yupper = 50.e0

    topography = Topography(topo_func=qinit)
    topography.x = np.linspace(xlower, xupper, nxpoints)
    topography.y = np.linspace(ylower, yupper, nypoints)
    topography.write(outfile, topo_type=1)

def topo(x,y):
    """
    Flat 
    """
    return -10*np.ones(np.shape(x))


def qinit(x,y):
    """
    Uniform X flow
    """
    return 13*np.ones(np.shape(x))

if __name__=='__main__':
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = Path()
    maketopo(path)
    makeqinit(path)
