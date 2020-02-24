
from __future__ import print_function
from pylab import *

def latlong(d,theta,phi,Rearth):
    """
    Take a point at distance d from the North pole with longitude theta
    and return longitude and latitude of point after rotating pole down 
    to latitude phi.
    Applying this function to a set of points on a circle of radius d 
    will result in points that are equi-distant (distance d on the earth) 
    from the point at latitude phi, longitude 0.  

    Used to construct a radially symmetric ocean on the earth and place 
    gauges, islands, etc.
    """

    # Convert phi to radians:
    theta = theta * pi/180.
    phi = phi * pi/180.

    alpha = pi/2. - d/Rearth   # latitude of original point
    # x,y,z coordinates of this point on the earth:
    x1 = cos(alpha)*cos(theta)
    y1 = cos(alpha)*sin(theta)
    z1 = sin(alpha)

    # rotate so centered at latitude phi:
    x2 = x1*sin(phi) + z1*cos(phi)
    y2 = y1
    z2 = -x1*cos(phi) + z1*sin(phi)

    # compute longitude xhat and latitude yhat:
    xhat = -arctan(y2/x2)
    yhat = arcsin(z2)

    # convert to degrees:
    xhat = xhat * 180./pi
    yhat = yhat * 180./pi

    return xhat,yhat

def plot_ocean_and_shelf(d1=1580e3, d2=1645e3, phi=40.):

    theta = linspace(0, 360., 200)
    d = d1*ones(theta.shape)
    xhat, yhat = latlong(d, theta, phi)

    clf()
    plot(xhat,yhat,'b')
    hold(True)

    d = d2*ones(theta.shape)
    xhat, yhat = latlong(d, theta, phi)

    plot(xhat,yhat,'r')
    legend(['Continental shelf', 'Shoreline'],loc='lower left')

    (xi1,yi1) = latlong(1600.e3,220.,40.)
    plot([xi1],[yi1],'ko')
    (xi2,yi2) = latlong(1600.e3,260.,40.)
    plot([xi2],[yi2],'ko')
    axis('scaled')
    xlim([-20, 20])
    ylim([15, 60])

#==========================================================
def gcdist(x1,y1,x2,y2,Rearth,units='degrees'):
    """
    Compute the great circle distance on the earth between points
    (x1,y1) and (x2,y2), where:
    x = longitude, y = latitude 

    Taken from Clawpack 4.x.
    """
    from numpy import pi,sin,cos,arccos,arcsin,sqrt
    if units=='degrees':
        # convert to radians:
        x1 = x1 * pi/180.
        y1 = y1 * pi/180.
        x2 = x2 * pi/180.
        y2 = y2 * pi/180.
    elif units != 'radians':
        raise Exception("unrecognized units")

    dx = x1 - x2
    dy = y1 - y2

    # angle subtended by two points, using Haversine formula:
    dsigma = 2. * arcsin(sqrt(sin(0.5*dy)**2 + cos(y1)*cos(y2)*sin(0.5*dx)**2))

    # alternative formula that may have more rounding error:
    #dsigma2 = arccos(sin(y1)*sin(y2)+ cos(y1)*cos(y2)*cos(dx))
    #print("max diff in dsigma: ", abs(dsigma-dsigma2).max())

    d = Rearth * dsigma
    return d


