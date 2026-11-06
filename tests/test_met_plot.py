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


# ---------------------------------------------------------------------------
#  Vorticity
# ---------------------------------------------------------------------------
def _cell_centers(num_cells, lower=0.0, upper=1.0):
    dx = (upper - lower) / num_cells
    x = lower + (np.arange(num_cells) + 0.5) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    return X, Y, dx


def test_vorticity_solid_body_rotation_exact():
    # u = -Omega y, v = Omega x has omega = 2 Omega; both the central and the
    # one-sided stencils are exact for linear fields, including at the edges.
    Omega = 0.7
    X, Y, dx = _cell_centers(10)
    omega = met_plot.extract_vorticity([dx, dx], -Omega * Y, Omega * X)
    np.testing.assert_allclose(omega, 2.0 * Omega, rtol=0, atol=1e-12)


def test_vorticity_second_order_interior_and_edges():
    # u = sin(x) sin(y), v = cos(x) cos(y) -> omega = -2 sin(x) cos(y)
    interior_err, edge_err = [], []
    sizes = (16, 32, 64, 128)
    for n in sizes:
        X, Y, dx = _cell_centers(n)
        omega = met_plot.extract_vorticity(
            [dx, dx], np.sin(X) * np.sin(Y), np.cos(X) * np.cos(Y))
        err = np.abs(omega + 2.0 * np.sin(X) * np.cos(Y))
        edge = np.zeros(err.shape, dtype=bool)
        edge[[0, -1], :] = True
        edge[:, [0, -1]] = True
        interior_err.append(err[~edge].max())
        edge_err.append(err[edge].max())

    for errors in (interior_err, edge_err):
        rates = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
        assert np.all(rates > 1.9), rates


@pytest.mark.parametrize("dry_cell, expected", [
    # interior: the cell and its four neighbors
    ((4, 4), {(4, 4), (3, 4), (5, 4), (4, 3), (4, 5)}),
    # next to the x = lower edge: the edge cell's one-sided stencil reaches it
    ((1, 4), {(0, 4), (1, 4), (2, 4), (1, 3), (1, 5)}),
])
def test_vorticity_dry_mask(dry_cell, expected):
    X, Y, dx = _cell_centers(8)
    dry = np.zeros(X.shape, dtype=bool)
    dry[dry_cell] = True
    omega = met_plot.extract_vorticity([dx, dx], -Y, X, dry=dry)

    assert set(zip(*np.nonzero(np.isnan(omega)))) == expected
    np.testing.assert_allclose(omega[~np.isnan(omega)], 2.0, atol=1e-12)


def test_water_vorticity_uses_patch_spacing():
    # Anisotropic spacing checks dx goes with x (axis 0) and dy with y
    from types import SimpleNamespace
    nx, ny, dx, dy, Omega = 6, 5, 0.1, 0.3, 1.5
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")
    h = 2.0 * np.ones(X.shape)
    h[0, 0] = 0.0
    q = np.stack([h, h * (-Omega * Y), h * (Omega * X), h])
    omega = met_plot.water_vorticity(SimpleNamespace(q=q, dx=dx, dy=dy))

    assert np.isnan(omega[0, 0]) and np.isnan(omega[1, 0])
    np.testing.assert_allclose(omega[2:, 1:], 2.0 * Omega, atol=1e-12)


@pytest.mark.parametrize("num_cells", [1, 2])
def test_vorticity_small_patches(num_cells):
    X, Y, dx = _cell_centers(num_cells)
    omega = met_plot.extract_vorticity([dx, dx], -Y, X)
    if num_cells == 1:
        assert np.all(np.isnan(omega))
    else:
        np.testing.assert_allclose(omega, 2.0, atol=1e-12)
