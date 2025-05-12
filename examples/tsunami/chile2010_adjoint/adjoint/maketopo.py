"""
    Download topo and dtopo files needed for this example.
    
    Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

from pathlib import Path
import os
import sys

import numpy as np

import clawpack.clawutil.data
import clawpack.geoclaw.topotools as topotools
from clawpack.geoclaw.util import haversine

try:
    CLAW = Path(os.environ['CLAW'])
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = CLAW / 'geoclaw' / 'scratch'


# Initial data for adjoint is Gaussian hump around this location:
# DART 32412 location:
DART_32412_location = (-86.392, -17.975)


def get_topo(path=None, makeplots=False):
    """
    Retrieve the topo file from the GeoClaw repository.
    """

    if not path:
        path = scratch_dir

    topo_fname = 'etopo10min120W60W60S0S.asc'
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/etopo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=path, 
            file_name=topo_fname, verbose=True)

    if makeplots:
        import matplotlib.pyplot as plt
        topo = topotools.Topography(path / topo_fname, topo_type=2)
        topo.plot()
        fname = topo_fname.with_suffix('.png')
        plt.savefig(fname)
        print("Created ", fname)


def makeqinit(path=None, center=None):
    """
        Create qinit data file
    """

    if not center:
        center = DART_32412_location
    if not path:
        path = Path()
    else:
        path = Path(path)

    nxpoints = 201
    nypoints = 201
    
    xlower = center[0] - 1.5
    xupper = center[0] + 1.5
    ylower = center[1] - 1.5
    yupper = center[1] + 1.5
    
    outfile = path / "hump.xyz"
    topotools.topo1writer(outfile, lambda x, y: qinit(x, y, center=center), 
                            xlower, xupper, ylower, yupper, nxpoints, nypoints)

def qinit(x, y, center=None):

    if not center:
        center = DART_32412_location

    # Gaussian using distance in meters:
    r = haversine(x, y, center[0], center[1])
    return np.exp(-(r / 20e3)**2)

if __name__=='__main__':
    get_topo(False)
    makeqinit()
