#!/usr/bin/env python

# import os
from pathlib import Path
import sys
import shutil
import pytest
from urllib.parse import urlparse
from urllib.error import URLError

import numpy as np

import clawpack.geoclaw.topotools as topotools
import clawpack.clawutil.data

# Set test data path
test_data_path = Path(__file__).parent / "test_data"

def topo_bowl(x,y):
    """Sample topo"""
    return 1000.0 * (x**2 + y**2 - 1.0)


def test_read_write_topo_bowl(tmp_path):
    """
    Test writing and reading topo files with small number of points
    Note that ordering should go from NW corner.
    """

    # Base topography
    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = np.linspace(-1.0, 3.0, 5)
    topo.y = np.linspace( 0.0, 3.0, 4)

    assert np.allclose(topo.x, np.array([-1.,  0.,  1.,  2., 3.])), \
           "Topography x values are incorrect."
    assert np.allclose(topo.X,
                          np.array([[-1.,  0.,  1.,  2.,  3.],
                                       [-1.,  0.,  1.,  2.,  3.],
                                       [-1.,  0.,  1.,  2.,  3.],
                                       [-1.,  0.,  1.,  2.,  3.]])), \
           "Topography X values are incorrect."
    assert np.allclose(topo.y, np.array([ 0.,  1.,  2.,  3.])), \
           "Topography y values are incorrect."
    assert np.allclose(topo.Y,
                          np.array([[ 0.,  0.,  0.,  0.,  0.],
                                       [ 1.,  1.,  1.,  1.,  1.],
                                       [ 2.,  2.,  2.,  2.,  2.],
                                       [ 3.,  3.,  3.,  3.,  3.]])), \
           "Topography Y values are incorrect."
    assert np.allclose(topo.Z,
                np.array([[     0.,  -1000.,      0.,   3000., 8000.],
                             [  1000.,      0.,   1000.,   4000.,   9000.],
                             [  4000.,   3000.,   4000.,   7000.,  12000.],
                             [  9000.,   8000.,   9000.,  12000.,  17000.]])), \
           "Topography Z values are incorrect."

    try:
        for topo_type in range(1, 4):
            path = Path(tmp_path) / f'bowl.tt{topo_type}'
            topo.write(path, topo_type=topo_type,Z_format="%22.15e")

            topo_in = topotools.Topography(path)
            assert np.allclose(topo.Z, topo_in.Z), \
                   "Differnece in written and read topography found."
    except Exception as e:
        # Copy output if something raised an exception
        local_path = Path() / "test_read_write_topo_bowl"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


def test_crop_topo_bowl():
    """
    Test cropping a topo file.
    """

    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = np.linspace(-1.0, 3.0, 5)
    topo.y = np.linspace( 0.0, 3.0, 4)

    # topo.Z should be created automatically when referenced below:
    assert np.allclose(topo.Z,
           np.array([[     0.,  -1000.,      0.,   3000., 8000.],
                        [  1000.,      0.,   1000.,   4000.,   9000.],
                        [  4000.,   3000.,   4000.,   7000.,  12000.],
                        [  9000.,   8000.,   9000.,  12000.,  17000.]])), \
           "Basic topography does not match test data."

    cropped_topo = topo.crop([0, 1, 0, 2])
    assert np.allclose(cropped_topo.x, np.array([0.0, 1.0])), \
           "Cropped topography y values do not match"
    assert np.allclose(cropped_topo.y, np.array([ 0.,  1.,  2.])), \
           "Cropped topography y values do not match."
    assert np.allclose(cropped_topo.Z,
                          np.array([[-1000.,     0.],
                                       [    0.,  1000.],
                                       [ 3000.,  4000.]])), \
           "Cropped topography Z values do not match."


def test_old_topotools(tmp_path):
    """
    Test against the old topotools from 5.1.0.
    Compare bowl.tt1 to bowl_old.tt1
    """

    old_topotools = pytest.importorskip("clawpack.geoclaw.old_topotools")

    nxpoints = 5
    nypoints = 4
    xlower = -1.0
    xupper = 3.0
    ylower = 0.0
    yupper = 3.0
    topo = topotools.Topography(topo_func=topo_bowl)
    topo.x = np.linspace(xlower, xupper, nxpoints)
    topo.y = np.linspace(ylower, yupper, nypoints)

    temp_path = tempfile.mkdtemp()
    try:
        file_path = os.path.join(temp_path, "bowl_old.tt1")
        old_topotools.topo1writer(file_path, topo_bowl, xlower, xupper, ylower, 
                                  yupper, nxpoints, nypoints)
        X, Y, Z = old_topotools.topofile2griddata(file_path, topotype=1)
        Y = np.flipud(Y)
        Z = np.flipud(Z)

        assert np.allclose(topo.X, X), "Difference in X grid."
        assert np.allclose(topo.Y, Y), "Difference in Y grid."
        assert np.allclose(topo.Z, Z), "Difference in Z grid."

    except Exception as e:
        # Copy output if something raised an exception
        local_path = Path() / "test_old_topotools"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


def topo_bowl_hill(x,y):
    """Add gaussian hill to topo_bowl"""
    return topo_bowl(x, y) + 1000.0 * np.exp(-100 * ((x - 0.7)**2 + (y - 0.8)**2))


def test_read_write_topo_bowl_hill(tmp_path):
    """
    Test writing and reading topo files.
    """
    try:
        topo = topotools.Topography(topo_func=topo_bowl_hill)
        topo.x = np.linspace(-1.5, 2.5, 101)
        topo.y = np.linspace(-1.0, 2.0, 76)

        for topo_type in range(1,4):
            file_path = Path(tmp_path) / f'bowl_hill.tt{topo_type}'
            topo.write(file_path, topo_type=topo_type,Z_format="%22.15e")
            topo_in = topotools.Topography(path=file_path, topo_type=topo_type)
            assert np.allclose(topo.Z, topo_in.Z),      \
                        (f"Written file of topo_type={topo_type} " + 
                          "does not equal read in file.")
    except Exception as e:
        # Copy output if something raised an exception
        local_path = Path() / "test_read_write_topo_bowl_hill"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


def test_netcdf(tmp_path):
    r"""Test Python NetCDF formatted topography reading"""

    netCDF4 = pytest.importorskip("netCDF4")

    try:
        url = ("https://raw.githubusercontent.com/rjleveque/geoclaw/" + 
               "5f675256c043e59e5065f9f3b5bdd41c2901702c/" +
               "src/python/geoclaw/tests/kahului_sample_1s.tt2")
        clawpack.clawutil.data.get_remote_file(url, output_dir=tmp_path, 
                                                    force=True)

        # Paths
        local_path = Path(tmp_path) / Path(urlparse(url).path).name
        nc_path = Path(tmp_path) / "test.nc"

        # Write out NetCDF version of file
        ascii_topo = topotools.Topography(path=local_path)
        ascii_topo.read()
        ascii_topo.write(nc_path, topo_type=4,Z_format="%22.15e")

        # Read back in NetCDF file
        nc_topo = topotools.Topography(path=nc_path)
        nc_topo.read()

        # Compare arrays - use tolerance based on 30 arcsecond accuracy
        assert np.allclose(ascii_topo.x, nc_topo.x), \
                    "Flat x-arrays did not match."
        assert np.allclose(ascii_topo.y, nc_topo.y), \
                    "Flat y-arrays did not match."
        assert np.allclose(ascii_topo.Z, nc_topo.Z), \
                    "Flat y-arrays did not match."

    except URLError as e:
        pytest.skip("Could not fetch remote file, skipping test.")
        raise e

    except Exception as e:
        # Copy output if something raised an exception not caught above
        local_path = Path() / "test_read_netcdf"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


def test_get_remote_file(tmp_path):
    """Test the ability to fetch a remote file from the web."""
    
    try:
        url = ("https://raw.githubusercontent.com/rjleveque/geoclaw/" + 
               "5f675256c043e59e5065f9f3b5bdd41c2901702c/" + 
               "src/python/geoclaw/tests/kahului_sample_1s.tt2")
        clawpack.clawutil.data.get_remote_file(url, output_dir=tmp_path,
                                                    force=True)

        local_path = Path(tmp_path) / Path(urlparse(url).path).name
        download_topo = topotools.Topography(path=local_path)

        test_path = Path(test_data_path) / Path(urlparse(url).path).name
        test_topo = topotools.Topography(path=test_path)

        assert np.allclose(download_topo.Z, test_topo.Z), \
               "Downloaded file does not match %s" % test_path

    except URLError:
        raise nose.SkipTest("Could not fetch remote file, skipping test.")

    except Exception as e:
        # Copy output if something raised an exception not caught above
        shutil.copy(local_path, Path() / "test_remote_file.tt2")
        raise e


def test_unstructured_topo(save=False, plot=False):
    """Test unstrucutred topography support"""

    scipy = pytest.importorskip("scipy")

    # Create random test data
    def test_topo(x, y):
        return x * (1 - x) * np.cos(4 * np.pi * x) * np.sin(4 * np.pi * y**2)**2

    fill_topo = topotools.Topography()
    fill_topo.x = np.linspace(0, 1, 100)
    fill_topo.y = np.linspace(0, 1, 200)
    fill_topo.Z = test_topo(fill_topo.X, fill_topo.Y)

    points = np.loadtxt(Path(test_data_path) / "unstructured_points.txt")
    values = test_topo(points[:,0], points[:,1])

    # Create topography object
    topo = topotools.Topography(unstructured=True)
    topo.x = points[:,0]
    topo.y = points[:,1]
    topo.z = values

    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(16,6))
        axes = fig.add_subplot(1, 3, 1)
        fill_topo.plot(axes=axes)
        axes.set_title("True Field")
        axes = fig.add_subplot(1, 3, 2)
        topo.plot(axes=axes, region_extent=[0, 1, 0, 1])
        axes.set_title("Unstructured Field")

    topo.interp_unstructured(fill_topo, extent=[0, 1, 0, 1], delta=(1e-2,1e-2))
    assert not topo.unstructured

    # Load (and save) test data and make the comparison
    path = Path(test_data_path) / "unstructured_test_data.tt3"
    if save:
        topo.write(path, Z_format="%22.15e")

    compare_data = topotools.Topography(path=path)

    assert np.allclose(compare_data.Z, topo.Z)

    if plot:
        axes = fig.add_subplot(1, 3, 3)
        topo.plot(axes=axes)
        axes.set_title("Interpolated Field")

        plt.show()


def test_topo_plot(tmp_path):
    """
    Create topo and write out, then read in again and plot.
    Note that center of bowl should be at (0,0).

    :TODO:
     - [] Add plot output test comparison
    """

    matplotlib = pytest.importorskip("matplotlib")
    # Use windowless frontend for image generation only
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        topo = topotools.Topography(topo_func=topo_bowl_hill)
        topo.x = np.linspace(-1.5, 2.5, 101)
        topo.y = np.linspace(-1.0, 2.0, 76)

        topo.plot()
        fname = Path(tmp_path) / "bowl_hill.png"
        plt.savefig(fname)
        topo2 = topo.crop([0.5, 1.5, 0., 2.])
        topo2.plot()
        plt.title("Cropped topography")
        fname = Path(tmp_path) / "bowl_hill_crop.png"
        plt.savefig(fname)

    except Exception as e:
        # Copy output if something raised an exception not caught above
        local_path = Path() / "test_topo_plot"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e


def test_plot_kahului(tmp_path):
    r"""
    Example illustrating reading in a topo file and plotting.
    Uses the test data kahului_sample_1s.tt2, created by cropping 
    the data file obtained from the NGDC site
        http://www.ngdc.noaa.gov/dem/squareCellGrid/download/604
    In addition to using the Topography.plot function, also 
    illustrate how to do a contour data of the data directly.

    :TODO:
     - [] Add plot output test comparison
    """

    matplotlib = pytest.importorskip("matplotlib")
    # Use windowless frontend for image generation only
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        path = Path(test_data_path) / 'kahului_sample_1s.tt2'
        K = topotools.Topography(path, topo_type=2)

        assert K.Z.shape == (46, 65), "K.Z is wrong shape"
        assert np.allclose(K.Z[:3,:3], np.array([[ 11.339,  11.339, 11.339],
                                                 [ 13.339,  11.339, 11.339],
                                                 [ 13.339,  11.339, 10.339]])),\
               "Topography K does not match"

        fig, ax = plt.subplots()
        K.plot(axes=ax)
        ax.set_title("Kahului Harbor at 1 second resolution")
        fig.savefig(Path(tmp_path) / "kahului_imshow.png")

        # Make a contour plot of topography / bathymetry:
        fig, ax = plt.subplots()
        ax.contour(K.X, K.Y, K.Z, np.linspace(-20,-2,10), colors='b',
                                                          linestyles='-')
        ax.contour(K.X, K.Y, K.Z, np.linspace(2,20,10), colors='g')
        ax.contour(K.X, K.Y, K.Z, [0.], colors='r')

        # fix aspect ratio based on latitude:
        mean_lat = 0.5 * (K.y.max() + K.y.min())
        ax.set_aspect(1.0 / np.cos(np.pi / 180.0 * mean_lat))

        # fix tick marks so readable:
        ax.ticklabel_format(style="plain", useOffset=False)
        plt.xticks(rotation=20)

        plt.title("2-meter contours of topo (green) and bathymetry (blue)",\
                  fontsize=12)
        plt.savefig(Path(tmp_path) / "kahului_contour.png")

    except Exception as e:
        # Copy output if something raised an exception not caught above
        local_path = Path() / "test_kahului_plot"
        shutil.rmtree(local_path, ignore_errors=True)
        shutil.copytree(tmp_path, local_path)
        raise e

# :TODO:
#  - [] Add CLI capability including saving output data and plotting
# if __name__ == "__main__":
#     if len(sys.argv) > 1:
#         if "plot" in sys.argv[1].lower():
#             plot_kahului()
#             plot_topo_bowl_hill()
#             test_unstructured_topo(save=False, plot=True)
#         elif bool(sys.argv[1]):
#             test_unstructured_topo(save=True)
#     else:
#         # Run tests one at a time
#         test_read_write_topo_bowl()
#         test_crop_topo_bowl()
#         test_against_old()
#         test_read_write_topo_bowl_hill()
#         test_get_remote_file()
#         test_unstructured_topo()
#         test_netcdf()

#         print("All tests passed.")
