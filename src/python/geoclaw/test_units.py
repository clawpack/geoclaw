#!/usr/bin/env python
# encoding: utf-8

import numpy as np

from clawpack.geoclaw.units import units, convert

def test_conversions(verbose=False):
    r"""Test unit conversions."""

    for (measurement_type, measurement_units) in units.items():
        value = np.pi
        units_list = list(units[measurement_type].keys())
        for i in range(len(units_list)):
            if verbose:
                print("%s (%s) -> (%s)" % (value, units_list[i - 1], 
                                                  units_list[i]))
            value = convert(value, units_list[i - 1], units_list[i])
        np.testing.assert_allclose([value], [np.pi],
                      err_msg="Measurement tyep %s failed." % measurement_type)

if __name__ == '__main__':
    test_conversions(verbose=True)
