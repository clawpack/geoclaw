"""
Create topo and dtopo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    dtopo_usgs100227.tt3              create using Okada model 
Prior to Clawpack 5.2.1, the fault parameters we specified in a .cfg file,
but now they are explicit below.
    
Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

from __future__ import absolute_import
from __future__ import print_function
import os

import clawpack.clawutil.data

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW environment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')

def get_topo(makeplots=False):
    """
    Retrieve the topo file from the GeoClaw repository.
    """
    from clawpack.geoclaw import topotools

    topo_fname = 'etopo1min130E210E0N60N.asc'
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/etopo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)

    topo_fname = 'hawaii_6s.txt'
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/hawaii/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)

    topo_fname = 'kahului_1s.txt'
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/hawaii/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)

    topo_fname = 'fujii.txydz'
    url = 'http://depts.washington.edu/clawpack/geoclaw/dtopo/tohoku/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)

    if makeplots:
        from matplotlib import pyplot as plt
        topo = topotools.Topography(os.path.join(scratch_dir,topo_fname), topo_type=2)
        topo.plot()
        fname = os.path.splitext(topo_fname)[0] + '.png'
        plt.savefig(fname)
        print("Created ",fname)



if __name__=='__main__':
    get_topo(False)
