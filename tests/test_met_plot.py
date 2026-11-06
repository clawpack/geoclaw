"""Tests for the water helper functions in clawpack.geoclaw.met.plot"""

import numpy as np
import pytest

pytest.importorskip("matplotlib")
from clawpack.geoclaw.met import plot as met_plot
from clawpack.geoclaw.multilayer import plot as multilayer_plot


def test_extract_eta_masks_dry_and_nan_depth():
    h = np.array([[1.0, 0.0], [np.nan, 2.0]])
    eta = np.array([[1.5, 0.5], [0.25, 3.0]])

    out = met_plot.extract_eta(h, eta.copy())

    # NaN depth was never masked: `h == np.nan` is always False
    assert np.isnan(out[1, 0])
    assert np.isnan(out[0, 1])
    assert out[0, 0] == 1.5 and out[1, 1] == 3.0


@pytest.mark.parametrize("module, dry_value", [(met_plot, 0.0),
                                               (multilayer_plot, np.nan)])
def test_extract_velocity_dry_and_nan_depth(module, dry_value):
    h = np.array([[2.0, 0.0], [np.nan, 4.0]])
    hu = np.array([[1.0, 1.0], [1.0, 2.0]])

    u = module.extract_velocity(h, hu)

    np.testing.assert_array_equal(u[[0, 1], [0, 1]], [0.5, 0.5])
    np.testing.assert_array_equal(u[[0, 1], [1, 0]], [dry_value, dry_value])
