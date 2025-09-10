"""
The function adjust_xy will center a point in a finite volume grid cell of
a given resolution on a given domain.

This can be used in particular to take adjust an approximate point
(x_desired, y_desired) where you want a computational gauge, to obtain
a point that lies exactly in the center of a grid cell at a given resolution.
This eliminates the need to interpolate between cell values in GeoClaw output,
which has issues for cells near the shoreline as described at
https://www.clawpack.org/nearshore_interp.html

Functions:

 - adjust: utility function to adjust in one dimension (x or y)
 - adjust_xy: the function to call to center a point in both x and y.
 - test: a simple example

"""

import numpy as np

def adjust(z_desired, z_edge, dz, verbose=False):
    """
    Given a desired location z (either x or y) for a gauge or other point
    of interest, adjust the location is it is offset (integer + 1/2)*dz
    from z_edge, an arbitrary edge location (e.g. an edge of the
    computational domain or any point offset integer*dz from the domain edge).
    This will put the new location in the center of a finite volume cell
    provided this point is in a patch with resolution dz.
    """
    # z represents either x or y
    i = np.round((z_desired-z_edge - 0.5*dz)/dz)
    z_centered = z_edge + (i+0.5)*dz
    if verbose:
        zshift = z_centered - z_desired
        zfrac = zshift / dz
        print("adjusted from %15.9f" % z_desired)
        print("           to %15.9f" % z_centered)
        print("   shifted by %15.9f = %.3f*dz" % (zshift,zfrac))
    return z_centered

def adjust_xy(x_desired, y_desired, x_edge, y_edge, dx, dy, verbose=False):
    """
    Given a desired location (or array of locations) in 2D, center
    so that the point(s) are at the centers of cells of size dx by dy,
    with cell edges that are integer multiples of dx,dy away from
    from x_edge,y_edge  (e.g. the edges of the computational domain or
    any other points offset by integer multiples of dx,dy from the edges.

    This will put the new location in the center of a finite volume cell
    provided this point is in a grid patch with resolution dx,dy.

    :Input:

    - x_desired, y_desired: single floats or arrays of floats, the desired
      locations
    - x_edge, y_edge (float) the edges of the computational domain
      (or any other points offset by integer multiples of dx,dy from the edges)
    - dx, dy (float): the grid resolution on which the point(s) should be
      centered

    :Output:

    - xc,yc: centered points that lie within dx/2, dy/2 of the desired
      location(s) and with (xc - x_edge)/dx and (yc - y_edge)/dy
      equal to an integer + 0.5,
    """

    # convert single values or lists/tuples to numpy arrays
    x_desired = np.array(x_desired, ndmin=1)
    y_desired = np.array(y_desired, ndmin=1)

    assert len(x_desired) == len(y_desired), \
            '*** lengths of x_desired, y_desired do not match'

    x_centered = []
    y_centered = []
    for i in range(len(x_desired)):
        x = adjust(x_desired[i], x_edge, dx, verbose=verbose)
        y = adjust(y_desired[i], y_edge, dy, verbose=verbose)
        x_centered.append(x)
        y_centered.append(y)

    if len(x_centered) == 1:
        return float(x_centered[0]), float(y_centered[0])
    else:
        return np.array(x_centered), np.array(y_centered)

def test():

    x_desired = [-122.01, -122.02]
    y_desired = [47.001, 47.002]

    # grid resolution on which to center point:
    dx = 1/(3*3600.)

    # lower left edge of computational domain:
    x_edge = -123.
    y_edge = 45.

    print('Desired x = ',x_desired)
    print('Desired y = ',y_desired)

    xc,yc = adjust_xy(x_desired,y_desired,x_edge,y_edge,dx,dx,verbose=True)

    print('Centered x = ',xc)
    print('Offsets in x in units of 1/3 arcsec: ', (xc-x_edge)*3*3600)
    print('Centered y = ',yc)
    print('Offsets in y in units of 1/3 arcsec: ', (yc-y_edge)*3*3600)

if __name__ == '__main__':
    test()
