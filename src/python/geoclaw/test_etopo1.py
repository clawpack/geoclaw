#!/usr/bin/env python

from pathlib import Path
import warnings

import numpy as np
import pytest

import clawpack.geoclaw.topotools as topotools

# Set test data path
test_data_path = Path(__file__).parent / "test_data"

def test_etopo1_topo(save=False, extent=[-125, -124, 48, 48.5]):
    """Test ETOPO1 support"""

    netCDF4 = pytest.importorskip("netCDF4")

    try:
        topo1 = topotools.read_netcdf('etopo1', extent=extent, verbose=True)
        topo10 = topotools.read_netcdf('etopo1', extent=extent, 
                                       coarsen=10, verbose=True)
    except (OSError, RuntimeError):
        warnings.warn('Could not read etopo1 data, check if thredds server up')
        pytest.skip("Reading etopo1 failed, skipping test")

    path = test_data_path / 'etopo1_10min.asc'
    if save:
        topo10.write(path, topo_type=3, Z_format='%.0f')

    topo10input = topotools.Topography()
    topo10input.read(path, topo_type=3)
    
    assert np.allclose(topo10.Z, topo10input.Z), \
           "topo10.Z does not agree with archived data"
    
    # if make_plot:
    #     import matplotlib.pyplot as plt
    #     plt.figure(figsize=(12,5))
    #     ax1 = plt.subplot(1,2,1)
    #     topo1.plot(axes=ax1)
    #     plt.title('1 minute etopo1 data')
    #     ax10 = plt.subplot(1,2,2)
    #     topo10.plot(axes=ax10)
    #     plt.title('10 minute etopo1 data')
    #     pname = 'etopo1_test_plot.png'
    #     plt.savefig(pname)
    #     print('Created %s' % pname)
    
def test_etopo1_xarray(extent=[-125, -124, 48, 48.5]):
    """Test xarray topography support"""

    xarray = pytest.importorskip("xarray")
        
    try:
        topo10, topo10_xarray = topotools.read_netcdf('etopo1', extent=extent, 
                                                      return_xarray=True,
                                                      coarsen=10, verbose=True)
    except (OSError, RuntimeError):
        warnings.warn('Could not read etopo1 data, check if thredds server up')
        pytest.skip("Reading etopo1 failed, skipping test")

    path = test_data_path / 'etopo1_10min.asc'
    topo10input = topotools.Topography()
    topo10input.read(path, topo_type=3)
    
    assert np.allclose(topo10_xarray['z'], topo10input.Z), \
           "topo10_xarray['z'] does not agree with archived data"
    

# :TODO:
#  - [] Add CLI capability including saving output data and plotting
# if __name__ == "__main__":
#     import sys
#     if len(sys.argv) > 1:
#         if "plot" in sys.argv[1].lower():
#             test_etopo1_topo(make_plot=True)
#         elif bool(sys.argv[1]):
#             test_etopo1_topo(save=True)
#     else:
#         # Run tests
#         test_etopo1_topo()
#         test_etopo1_xarray()
#         print("All tests passed.")

