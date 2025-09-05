#!/usr/bin/env python

from pathlib import Path
import pytest

import numpy as np

import clawpack.geoclaw.surge.storm
import clawpack.geoclaw.util as util

file_format_map = {1: ["ascii", 1, "nws12", "owi"], 
                   2: ["netcdf", 2, "nws13"]}

def create_netcdf_storm_file(path):
    """"""
    xr = pytest.importorskip("xarray")
    size = (3, 5, 3)
    coords = [('longitude', np.linspace(-1, 1, size[0])), 
              ('latitude', np.linspace(-2, 2, size[1])), 
              ('valid_time', np.linspace(0, 2, size[2]))]
    wind_x = xr.DataArray(np.random.rand(*size),  coords=coords)
    wind_y = xr.DataArray(np.random.rand(*size), coords=coords)
    P = xr.DataArray(np.random.rand(*size), coords=coords)
    ds = xr.Dataset({'u': wind_x, 'v': wind_y, 'pressure': P,})
    ds.to_netcdf(path)


@pytest.mark.parametrize("file_format", ['ascii', 'netcdf'])
def test_data_storms(file_format, tmp_path):
    """Test of Storm OWI formatted I/O"""

    storm_path = Path(tmp_path) / "test.storm"
    storm = clawpack.geoclaw.surge.storm.Storm()
    storm.time_offset = np.datetime64("2012-08-29")
    storm.file_format = file_format
    storm.window_type = 1
    storm.ramp_width = 3
    if file_format == 'ascii':
        storm.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                            Path("storm_2.PRE"), Path("storm_2.WIN")]
        storm.write(storm_path, file_format="data")
        read_storm = clawpack.geoclaw.surge.storm.Storm(storm_path, 
                                                        file_format="data")
    else:
        storm.file_paths = [tmp_path / "storm.nc"]
        create_netcdf_storm_file(storm.file_paths[0])
        storm.write(storm_path, file_format="data", 
                                dim_mapping={"t": "valid_time"})
        read_storm = clawpack.geoclaw.surge.storm.Storm(storm_path, 
                                                        file_format="data")

    assert (storm.time_offset == read_storm.time_offset)
    assert (storm.file_format in
            file_format_map[read_storm.file_format])
    assert (storm.window_type == read_storm.window_type)
    assert (storm.ramp_width == read_storm.ramp_width)
    for (i, path) in enumerate(storm.file_paths):
        assert read_storm.file_paths[i] == path


def test_netcdf_var_mapping(tmp_path):
    """Test NetCDF file reading"""

    storm_data_file = Path(tmp_path) / "storm.nc"
    create_netcdf_storm_file(storm_data_file)

    # Test name finding
    _dim_mapping = util.get_netcdf_names(storm_data_file, 
                                         lookup_type='dim',
                                         verbose=True,
                                         user_mapping={'t': 'valid_time'})
    _var_mapping = util.get_netcdf_names(storm_data_file, 
                                         lookup_type='var',
                                         verbose=True)
    assert _dim_mapping == {'x': 'longitude', 
                            'y': 'latitude', 
                            't': 'valid_time'}
    assert _var_mapping == {'wind_u': 'u', 
                            'wind_v': 'v', 
                            'pressure': 'pressure'}
