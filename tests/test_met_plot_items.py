"""Tests for the plot item helpers (add_*) in clawpack.geoclaw.met.plot"""

import numpy as np
import pytest

pytest.importorskip("matplotlib")
from clawpack.visclaw.data import ClawPlotData
from clawpack.geoclaw.met import plot as met_plot


@pytest.fixture
def plotaxes():
    return ClawPlotData().new_plotfigure(name="test").new_plotaxes()


def _item(plotaxes, name):
    return plotaxes.plotitem_dict[name]


def test_surface_elevation_contourf_bounds(plotaxes):
    # Used to raise NameError (numpy.linspace without numpy imported)
    met_plot.add_surface_elevation(plotaxes, plot_type="contourf",
                                   bounds=[-1.0, 1.0])
    item = _item(plotaxes, "surface")
    np.testing.assert_allclose(item.contour_levels, np.linspace(-1, 1, 11))
    assert (item.fill_cmin, item.fill_cmax) == (-1.0, 1.0)
    assert item.fill_cmap is met_plot.surface_cmap


def test_surface_elevation_contourf_defaults(plotaxes):
    # Used to raise TypeError iterating contour_levels = None
    met_plot.add_surface_elevation(plotaxes, plot_type="contourf")
    assert _item(plotaxes, "surface").contour_levels is None


def test_speed_contourf_bounds(plotaxes):
    # Used to raise NameError, and then set fill_cmap = bounds[1]
    met_plot.add_speed(plotaxes, plot_type="contourf", bounds=[0.0, 2.0])
    item = _item(plotaxes, "speed")
    np.testing.assert_allclose(item.contour_levels, np.linspace(0, 2, 11))
    assert (item.fill_cmin, item.fill_cmax) == (0.0, 2.0)
    assert item.fill_cmap.name == "PuBu"


def test_wind_contour_levels(plotaxes):
    # contour_min was misspelled, and visclaw ignores contour_min without
    # contour_max, so the refinement thresholds were never used as levels
    met_plot.add_wind(plotaxes, plot_type="contour")
    item = _item(plotaxes, "wind")
    assert item.contour_levels == list(met_plot.surge_data.wind_refine)
    assert not hasattr(item, "countour_min")
