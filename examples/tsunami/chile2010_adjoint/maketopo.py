"""
Create topo and dtopo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    dtopo_usgs100227.tt3              create using Okada model 
Prior to Clawpack 5.2.1, the fault parameters we specified in a .cfg file,
but now they are explicit below.
    
Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

from pathlib import Path
import os

import numpy as np

import clawpack.clawutil.data
import clawpack.geoclaw.topotools as topotools
import clawpack.geoclaw.dtopotools as dtopotools

try:
    CLAW = Path(os.environ['CLAW'])
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = CLAW / 'geoclaw' / 'scratch'

def get_topo(path=None, verbose=False, makeplots=False):
    """
    Retrieve the topo file from the GeoClaw repository.
    """

    if not path:
        path = scratch_dir

    topo_fname = 'etopo10min120W60W60S0S.asc'
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/etopo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=path, 
            file_name=topo_fname, verbose=verbose)

    if makeplots:
        import matplotlib.pyplot as plt
        topo = topotools.Topography(path / topo_fname, topo_type=2)
        topo.plot()
        fname = topo_fname.with_suffix('.png')
        plt.savefig(fname)
        if verbose:
            print("Created ", fname)

    
def make_dtopo(path=None, verbose=False, makeplots=False):
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified below.
    """
    
    if not path:
        path = scratch_dir

    dtopo_fname = path / "dtopo_usgs100227.tt3"

    # Specify subfault parameters for this simple fault model consisting
    # of a single subfault:

    usgs_subfault = dtopotools.SubFault()
    usgs_subfault.strike = 16.
    usgs_subfault.length = 450.e3
    usgs_subfault.width = 100.e3
    usgs_subfault.depth = 35.e3
    usgs_subfault.slip = 15.
    usgs_subfault.rake = 104.
    usgs_subfault.dip = 14.
    usgs_subfault.longitude = -72.668
    usgs_subfault.latitude = -35.826
    usgs_subfault.coordinate_specification = "top center"

    fault = dtopotools.Fault()
    fault.subfaults = [usgs_subfault]

    if verbose:
        print("Mw = ",fault.Mw())

    if dtopo_fname.exists():
        if verbose:
            print("Not regenerating dtopo file (already exists): " + 
                  f"{dtopo_fname}")
    else:
        if verbose:
            print("Using Okada model to create dtopo file")

        x = np.linspace(-77, -67, 100)
        y = np.linspace(-40, -30, 100)
        times = [1.]

        fault.create_dtopography(x,y,times)
        dtopo = fault.dtopo
        dtopo.write(dtopo_fname, dtopo_type=3)

    if makeplots:
        import matplotlib.pyplot as plt
        if fault.dtopo is None:
            # read in the pre-existing file:
            print("Reading in dtopo file...")
            dtopo = dtopotools.DTopography()
            dtopo.read(dtopo_fname, dtopo_type=3)
            x = dtopo.x
            y = dtopo.y
        plt.figure(figsize=(12,7))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
        fault.plot_subfaults(axes=ax1,slip_color=True)
        ax1.set_xlim(x.min(),x.max())
        ax1.set_ylim(y.min(),y.max())
        dtopo.plot_dZ_colors(1.,axes=ax2)
        fname = dtopo_fname.with_suffix('.png')
        plt.savefig(fname)
        if verbose:
            print("Created ",fname)


if __name__=='__main__':
    get_topo()
    make_dtopo()
