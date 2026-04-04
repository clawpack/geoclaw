#!/usr/bin/env python
# encoding: utf-8

import pytest
import numpy as np

from clawpack.geoclaw.units import units, convert

@pytest.mark.python
@pytest.mark.parametrize("measurement_type, measurement_units", units.items())
def test_conversions(measurement_type, measurement_units):
    """Test unit conversions."""

    value = np.pi
    units_list = list(measurement_units.keys())
    for i in range(len(units_list)):
        value = convert(value, units_list[i - 1], units_list[i])
    np.testing.assert_allclose(
        value,
        np.pi,
        err_msg=f"Measurement type {measurement_type} failed.",
        )

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
