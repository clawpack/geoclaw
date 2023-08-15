#!/usr/bin/env python
# encoding: utf-8

"""Tests for reading and writing GeoClaw data files"""

import tempfile
import shutil
import os

import numpy

import clawpack.clawutil.test as test
import clawpack.geoclaw.data


def test_read_friction_data():
    r"""Test readinga and writing of FrictionData files"""

    # Test Data
    test_regions = []
    test_regions.append([(-99.0, -70.0), (8.0, 32.0), 
                       [numpy.infty, 0.0, -numpy.infty], 
                       [0.030, 0.022]])
    test_regions.append([(-98, 25.25), (-90, 30), 
                       [numpy.infty, -10.0, -200.0, -numpy.infty],
                       [0.030, 0.012, 0.022]])

    # Create temp directory
    temp_path = tempfile.mkdtemp()

    try:
        data_file = os.path.join(temp_path, "friction.data")
        
        # Create and write out data object
        friction_data = clawpack.geoclaw.data.FrictionData()
        friction_data.variable_friction = True
        for i in range(2):
            friction_data.friction_regions.append(test_regions[i])
        friction_data.write(data_file)

        # Read data object
        read_friction_data = clawpack.geoclaw.data.FrictionData()
        read_friction_data.read(data_file)

        # Tests
        for (i, region) in enumerate(test_regions):
            for j in range(4):
                assert numpy.allclose(region[j], 
                                      read_friction_data.friction_regions[i][j])

    except AssertionError as e:
        # If the assertion failed then copy the contents of the directory
        shutil.copytree(temp_path, os.path.join(os.getcwd(), 
                                                "test_read_friction_data"))
        raise e

    finally:
        shutil.rmtree(temp_path)


if __name__ == "__main__":
    test_read_friction_data()
    print("All tests passed.")