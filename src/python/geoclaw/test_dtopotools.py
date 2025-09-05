#!/usr/bin/env python

from pathlib import Path
import sys
import shutil
# import inspect
# import time
import pytest

import numpy as np

import clawpack.geoclaw.dtopotools as dtopotools

# Set test data path
test_data_path = Path(__file__).parent / "test_data"

def test_read_csv_make_dtopo(save=False):
    r"""Test reading and making of a CSV subfault speficied dtopo."""
    
    subfault_path = test_data_path / 'alaska1964.csv'
    input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
                   "mu":"dyne/cm^2"}
    fault = dtopotools.CSVFault()
    fault.read(subfault_path, input_units=input_units, 
                              coordinate_specification="noaa sift")

    assert abs(fault.Mw() - 8.53336) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 203
    xupper = 214.  # approximate - adjusted below
    ylower = 55
    yupper = 60.  # approximate - adjusted below

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1. / points_per_degree
    mx = int((xupper - xlower) / dx + 1)
    xupper = xlower + (mx - 1) * dx
    my = int((yupper - ylower) / dx + 1)
    yupper = ylower + (my - 1) * dx

    x = np.linspace(xlower, xupper, mx)
    y = np.linspace(ylower, yupper, my)

    dtopo = fault.create_dtopography(x, y, times=[1.])

    path = test_data_path /  "alaska1964_test_data.tt3"
    if save:
        dtopo.write(path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=path)
    compare_data.read(dtopo_type=3)

    assert dtopo.dZ.shape == compare_data.dZ.shape,     \
        (f"dtopo.dZ.shape is {dtopo.dZ.shape}, " + 
         f"should be {compare_data.dZ.shape}")

    assert np.allclose(compare_data.dZ, dtopo.dZ)


def test_read_ucsb_make_dtopo(save=False):
    r"""Test reading and making of a UCSB subfault speficied dtopo."""

    subfault_path = test_data_path / 'tohoku_ucsb.txt'
    fault = dtopotools.UCSBFault()
    fault.read(subfault_path)

    assert abs(fault.Mw() - 9.13957973) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 140.
    xupper = 146.
    ylower = 35.
    yupper = 41.

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1. / points_per_degree
    mx = int((xupper - xlower) / dx + 1)
    xupper = xlower + (mx - 1) * dx
    my = int((yupper - ylower) / dx + 1)
    yupper = ylower + (my - 1) * dx

    x = np.linspace(xlower, xupper, mx)
    y = np.linspace(ylower, yupper, my)

    tmax = 0.
    for s in fault.subfaults:
        tmax = max(tmax, s.rupture_time + s.rise_time)

    fault.rupture_type = 'kinematic'
    times = np.linspace(0, tmax, 10)
    dtopo = fault.create_dtopography(x, y, times)

    path = test_data_path / "tohoku_test_data.tt3"
    if save:
        dtopo.write(path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=path)
    compare_data.read(path=path, dtopo_type=3)

    assert dtopo.dZ.shape == compare_data.dZ.shape, \
        "dtopo.dZ.shape is %s, should be %s" \
        % (dtopo.dZ.shape, compare_data.dZ.shape)

    assert np.allclose(compare_data.dZ, dtopo.dZ)


def test_read_sift_make_dtopo(save=False):
    r"""Test reading and making of a SIFT subfault speficied dtopo"""

    sift_slip = {'acsza1':2, 'acszb1':3}
    fault = dtopotools.SiftFault(sift_slip)

    assert abs(fault.Mw() - 7.966666666667) < 1e-4, "*** Mw is wrong: %g" % fault.Mw()

    xlower = 162.
    xupper = 168.
    ylower = 53.
    yupper = 59.

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = np.linspace(xlower,xupper,mx)
    y = np.linspace(ylower,yupper,my)

    times = [1.]
    dtopo = fault.create_dtopography(x,y,times)

    path = test_data_path / "sift_test_data.tt3"
    if save:
        dtopo.write(path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=path)
    compare_data.read(dtopo_type=3)

    assert dtopo.dZ.shape == compare_data.dZ.shape, \
        "dtopo.dZ.shape is %s, should be %s" \
        % (dtopo.dZ.shape, compare_data.dZ.shape)

    assert np.allclose(compare_data.dZ, dtopo.dZ)


def test_SubdividedPlaneFault_make_dtopo(save=False):
    r"""Subdivided Plane Fault dtopo tests"""

    # get a unit source fault plane as starting point:
    sift_slip = {'acsza1':1.}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]
    # Mo = fault_plane.Mo()
    # print "original Mo = ",Mo

    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    # print "new Mo = ",fault2.Mo()
    #fault2.plot_subfaults(slip_color=True)

    assert abs(fault2.Mw() - 7.500686667) < 1e-4, \
           "*** Mw is wrong: %g" % fault.Mw()

    xlower = 162.
    xupper = 168.
    ylower = 53.
    yupper = 59.

    # dtopo parameters:
    points_per_degree = 4  # 15 minute resolution
    dx = 1./points_per_degree
    mx = int((xupper - xlower)/dx + 1)
    xupper = xlower + (mx-1)*dx
    my = int((yupper - ylower)/dx + 1)
    yupper = ylower + (my-1)*dx

    x = np.linspace(xlower,xupper,mx)
    y = np.linspace(ylower,yupper,my)

    times = [1.]
    dtopo = fault2.create_dtopography(x,y,times)

    path = test_data_path / "SubdividedFaultPlane_test_data.tt3"
    if save:
        dtopo.write(path, dtopo_type=3)
    compare_data = dtopotools.DTopography(path=path)
    compare_data.read(path=path, dtopo_type=3)

    assert dtopo.dZ.shape == compare_data.dZ.shape, \
        "dtopo.dZ.shape is %s, should be %s" \
        % (dtopo.dZ.shape, compare_data.dZ.shape)

    assert np.allclose(compare_data.dZ, dtopo.dZ)


@pytest.mark.skip(reason="Strangely not quite working in the dZ comparison")
def test_dtopo_io(tmp_path):
    r"""Test IO of dtopography class"""

    test_dtopo = dtopotools.DTopography(
                        path=test_data_path / "alaska1964_test_data.tt3")

    # temp_path = tempfile.mkdtemp()
    try:
        dtopo_paths = [Path(tmp_path) / 'alaska1964.tt1',
                       Path(tmp_path) / 'alaska1964.tt3']

        for path in dtopo_paths:
            test_dtopo.write(path)
            dtopo = dtopotools.DTopography(path=path)

            assert test_dtopo.dZ.shape == dtopo.dZ.shape,   \
                   f"Shape of dZ not equal for topo_type = {dtopo.topo_type}."

            assert np.allclose(test_dtopo.dZ, dtopo.dZ),    \
                   f"dZ not equal for {path}"

    except Exception as e:
        # Copy output if something raised an exception
        local_path = Path() / "test_dtopo_io"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e

# Test script for plotting differences
# #!/usr/bin/env python

# import os
# from pathlib import Path
# import matplotlib.pyplot as plt
# import clawpack.geoclaw.dtopotools as dtopotools
# test_data_path = (Path(os.environ["CLAW"]) / "geoclaw" / "src" / "python" / "geoclaw" / "test_data").resolve()

# def test_dtopo_io(tmp_path):
#     test_dtopo = dtopotools.DTopography(
#                         path=test_data_path / "alaska1964_test_data.tt3")
#     dtopo_paths = [Path(tmp_path) / 'alaska1964.tt1',
#                    Path(tmp_path) / 'alaska1964.tt3']

#     for path in dtopo_paths:
#         test_dtopo.write(path)
#         dtopo = dtopotools.DTopography(path=path)
#         print(test_dtopo.dZ.shape)
#         print(dtopo.dZ.shape)

#         fig, ax = plt.subplots()
#         plot = ax.pcolor(test_dtopo.dZ[0, :, :] - dtopo.dZ[0, :, :])
#         fig.colorbar(plot)
        
#     plt.show()

# if __name__ == '__main__':
#     test_dtopo_io(Path())


def test_geometry():
    r"""Test subfault geometry calculation."""

    old_dtopotools = pytest.importorskip("clawpack.geoclaw.old_dtopotools")

    subfault_path = test_data_path / 'alaska1964.csv'
    input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
                   "mu":"dyne/cm^2"}
    specifications = ['top center', 'centroid', 'bottom center', 'noaa sift']

    for specification in specifications:
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                                  coordinate_specification=specification)

        # Subfault 10 is chosen at random, maybe do all?
        subfault = fault.subfaults[10]
        geometry = old_dtopotools.set_geometry(subfault)

        coord_tests = {"top center":{'test':[geometry['x_top'], 
                                             geometry['y_top'], 
                                             geometry['depth_top']], 
                                 'computed':subfault.centers[0]},
                       "centroid":{'test':[geometry['x_centroid'],
                                           geometry['y_centroid']],
                               'computed':subfault.centers[1][:2]},
                       "bottom center":{"test":[geometry['x_bottom'],
                                                geometry['y_bottom'],
                                                geometry['depth_bottom']],
                                    "computed":subfault.centers[2]},
                       "Corner A":{"test":[geometry["x_corners"][2],
                                           geometry["y_corners"][2]],
                               "computed":subfault.corners[0][:2]},
                       "Corner B":{"test":[geometry["x_corners"][3],
                                           geometry["y_corners"][3]],
                               "computed":subfault.corners[1][:2]},
                       "Corner C":{"test":[geometry["x_corners"][0],
                                           geometry["y_corners"][0]],
                               "computed":subfault.corners[2][:2]},
                       "Corner D":{"test":[geometry["x_corners"][1],
                                           geometry["y_corners"][1]],
                               "computed":subfault.corners[3][:2]}

                    }

        for (values, coord_test) in coord_tests.items():
            assert np.allclose(coord_test['test'], coord_test['computed']), \
                   "Specification = %s, coords= %s:\n%s !=\n%s" % (
                                                         specification, 
                                                         values, 
                                                         coord_test['test'], 
                                                         coord_test['computed'])

def test_vs_old_dtopo(tmp_path):
    r"""Test new dtopotools with old version from 5.2"""

    old_dtopotools = pytest.importorskip("clawpack.geoclaw.old_dtopotools")

    try:
        subfault_path = test_data_path / 'alaska1964.csv'
        old_dtopo_file = Path(tmp_path) / 'old_alaska.tt1'
        old_dtopotools.builddynamicdeffile(subfault_path, subfault_path, 
                                           old_dtopo_file)
    
        input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m", 
                       "mu":"dyne/cm^2"}
        fault = dtopotools.CSVFault()
        fault.read(subfault_path, input_units=input_units, 
                                  coordinate_specification="noaa sift")
        new_dtopo = fault.create_dtopography()

        X, Y, dZ = old_dtopotools.read_dtopo_old(old_dtopo_file, 
                                                 deftype='dynamic', 
                                                 only_last=False)

        assert np.allclose(X, new_dtopo.X), \
               "X values do not agree between old and new dtopotools."

        assert np.allclose(Y, new_dtopo.Y), \
               "Y values do not agree between old and new dtopotools."

        assert np.allclose(dZ, new_dtopo.dZ), \
               "dZ values do not agree between old and new dtopotools."


    except Exception as e:
        # Copy output if something raised an exception
        local_path = Path() / "test_old_dtopo"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


def test_dynamic_tohoku(tmp_path):
    r"""Test dynamic faulting via a Tohoku example"""

    fault = dtopotools.UCSBFault()
    fault.read(test_data_path / 'tohoku_ucsb.txt')
    fault.rupture_type = 'dynamic'
    
    # Seafloor test
    xlower = 140.
    xupper = 146.
    ylower = 35.
    yupper = 41.
    xylim = [xlower, xupper, ylower, yupper]
    mx = int((xupper - xlower) * 15 + 1)
    my = int((yupper - ylower) * 15 + 1)
    x = np.linspace(xlower, xupper, mx)
    y = np.linspace(ylower, yupper, my)
    
    tmax = 0.
    for s in fault.subfaults:
        tmax = max(tmax, s.rupture_time + s.rise_time + s.rise_time_ending)
    times = np.linspace(0,tmax,10)
    dtopo = fault.create_dtopography(x,y,times,verbose=True)

    dz_final = dtopo.dZ[-1]
    dz_max = dz_final.max()

    # Test writing
    try:
        dtopo.write(Path(tmp_path) / 'tohoku_ucsb_dynamic.tt3', 3)
    except Exception as e:
        # Copy output if something raised an exception
        local_path = Path() / "test_dynamic_tohoku"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


# :TODO: Enable plotting tests
# def test_dynamic_plotting(tmp_path):
#     """Test plotting of dynamic fault"""

#     shoreline_xy = np.load(test_data_path / 'tohoku_shoreline_1min.npy')

#     fault = dtopotools.UCSBFault()
#     fault.read(test_data_path / 'tohoku_ucsb.txt')
#     fault.rupture_type = 'dynamic'

#     fig, ax = plt.subplots()
#     fault.plot_subfaults(axes=ax, slip_color=True)

#     # Incorporate this function in dtopotools to replace animate_dz_colors?
#     def plot_subfaults_dz(t, fig=None):
#         if fig is None:
#             fig = plt.figure(figsize=(12,5))
#         else:
#             fig.clf()
#         ax1 = fig.add_subplot(121)
#         ax2 = fig.add_subplot(122)
#         fault.plot_subfaults(axes=ax1, slip_color=True, slip_time=t, xylim=xylim)
#         dtopo.plot_dz_colors(axes=ax2, t=t, cmax_dz=dz_max)
#         ax1.plot(shoreline_xy[:,0],shoreline_xy[:,1],'g')
#         ax2.plot(shoreline_xy[:,0],shoreline_xy[:,1],'g')
#         plt.axis(xylim)
#         fig.show()
#         return fig

#     fig, ax = plt.subplots()
#     dtopo.plot_dz_colors(axes=ax, t=tmax)
#     plt.show()
            
#     fig = plt.figure(figsize=(12,5))

#     for t in list(np.linspace(0,150,16)) + [170,200]:
#         plot_subfaults_dz(t,fig)
#         plt.draw()
#         plt.show()
#         time.sleep(1)


def test_subdivided_plane_fault():
    r"""Test SubdividedPlaneFault class"""

    # get a unit source fault plane as starting point:
    sift_slip = {'acsza1':1.}
    fault = dtopotools.SiftFault(sift_slip)
    fault_plane = fault.subfaults[0]
    Mo = fault_plane.Mo()
    # if verbose:
    #     print("original Mo = ",Mo)

    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
    # if verbose:
    #     print("new Mo = ",fault2.Mo())
    # if plot:
    #     import matplotlib.pyplot as plt
    #     fault2.plot_subfaults(slip_color=True)
    #     plt.show()

    slip_function = lambda xi,eta: xi*(1-xi)*eta
    fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3, 
            slip_function=slip_function, Mo=Mo)
    # if verbose:
    #     print("new Mo = ",fault2.Mo())
    # if plot:
    #     fault2.plot_subfaults(slip_color=True)
    #     plt.show()

    fault2.subdivide(nstrike=20, ndip = 10, slip_function=slip_function, Mo=Mo)
    # if verbose:
    #     print("with finer resolution, Mo = ",fault2.Mo())
    # if plot:
    #     fault2.plot_subfaults(slip_color=True)
    #     plt.show()


# :TODO:
#  - [] Add CLI capability including saving output data and plotting
# if __name__ == "__main__":

#     save = False  # default 

#     if len(sys.argv) > 1:
#         if "plot" in sys.argv[1].lower():
#             pass
#         elif bool(sys.argv[1]):
#             save = True
#     try:
#         test_read_csv_make_dtopo(save=save)
#         test_read_ucsb_make_dtopo(save=save)
#         test_read_sift_make_dtopo(save=save)
#         test_SubdividedPlaneFault_make_dtopo(save=save)
#         test_dtopo_io()
#         test_geometry()
#         test_vs_old_dtopo()
#     except nose.SkipTest as e:
#         print(e.message)
