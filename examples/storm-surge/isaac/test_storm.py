#!/usr/bin/env python

import os
import datetime

import clawpack.geoclaw.surge.storm

test_paths = [os.path.join(os.getcwd(), "storm_1.PRE"),
              os.path.join(os.getcwd(), "storm_1.WIN"),
              os.path.join(os.getcwd(), "storm_2.PRE"),
              os.path.join(os.getcwd(), "storm_2.WIN")
             ]

for file_format in ['ascii', 'NWS12']:
    # Writing
    storm = clawpack.geoclaw.surge.storm.Storm()
    storm.time_offset = datetime.datetime(2012, 8, 29, 0)
    storm.data_file_format = file_format
    for path in test_paths:
        storm.file_paths.append(path)
    storm.write("test.storm", file_format="OWI")

    # Reading
    storm = clawpack.geoclaw.surge.storm.Storm("test.storm", file_format="OWI")
    assert storm.file_paths == test_paths

for file_format in ['NetCDF', 'NWS13']:
    # Writing
    storm = clawpack.geoclaw.surge.storm.Storm()
    storm.time_offset = datetime.datetime(2012, 8, 29, 0)
    storm.data_file_format = file_format
    storm.file_paths.append(test_paths[0])
    storm.write("test.storm", file_format="OWI")

    # Reading
    storm = clawpack.geoclaw.surge.storm.Storm("test.storm", file_format="OWI")
    assert storm.file_paths[0] == test_paths[0]