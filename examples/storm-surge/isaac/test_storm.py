#!/usr/bin/env python

from pathlib import Path
import pytest

import numpy as np

import clawpack.geoclaw.surge.storm

data_file_format_map = {1: ["ascii", 1, "nws12"], 
                        2: ["netcdf", 2, "nws13"]}

@pytest.mark.skip("Unimplemented test.")
def test_geoclaw_storm():
    """Test of Storm geoclaw formatted I/O"""
    assert False

@pytest.mark.parametrize("data_file_format", 
                         ['ascii', 'netcdf', 1, 2, 'nws12', 'nws13'])
def test_OWI_storm(data_file_format):
    """Test of Storm OWI formatted I/O"""
    storm = clawpack.geoclaw.surge.storm.Storm()
    storm.time_offset = np.datetime64("2012-08-29")
    storm.data_file_format = data_file_format
    if data_file_format in ['ascii', 1, 'nws12']:
        storm.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                            Path("storm_2.PRE"), Path("storm_2.WIN")]
    else:
        storm.file_paths = [Path("storm.nc")]
    storm.write(Path("test.storm"), file_format="OWI")
    read_storm = clawpack.geoclaw.surge.storm.Storm(Path("test.storm"), 
                                                    file_format="OWI")
    assert (storm.time_offset == read_storm.time_offset)
    assert (storm.data_file_format in
            data_file_format_map[read_storm.data_file_format])
    assert (storm.file_paths == read_storm.file_paths)

if __name__ == '__main__':
    test_OWI_storm()
