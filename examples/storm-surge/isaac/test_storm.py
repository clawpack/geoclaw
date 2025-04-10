#!/usr/bin/env python

from pathlib import Path
import numpy as np
import clawpack.geoclaw.surge.storm

def test_OWI_ascii_storm():
    """Test I/O for OWI ASCII/NWS12 format"""

    test_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                  Path("storm_2.PRE"), Path("storm_2.WIN")]

    storm = clawpack.geoclaw.surge.storm.Storm()
    storm.time_offset = np.datetime64("2012-08-29")
    storm.data_file_format = "ascii"
    for path in test_paths:
        storm.file_paths.append(path)
    storm.write("test.storm", file_format="OWI")

    storm = clawpack.geoclaw.surge.storm.Storm("test.storm", file_format="OWI")
    assert storm.file_paths == test_paths

def test_OWI_netcdf_storm():
    """Test I/O for OWI ASCII/NWS12 format"""

    test_paths = [Path("storm_1.nc")]

    storm = clawpack.geoclaw.surge.storm.Storm()
    storm.time_offset = np.datetime64("2012-08-29")
    storm.data_file_format = "netcdf"
    for path in test_paths:
        storm.file_paths.append(path)
    storm.write("test.storm", file_format="OWI")

    storm = clawpack.geoclaw.surge.storm.Storm("test.storm", file_format="OWI")
    assert storm.file_paths == test_paths
