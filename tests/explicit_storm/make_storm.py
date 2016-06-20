# encoding: utf-8
"""
Generate storm and bathymetry files
 for this very simple test case

This code generates a simple stationary storm field 
 based on the Holland model. It then writes the values
 to several data files meant as input in an 
 explicit storm module.

This is not intended to represent
 a realistic storm.
"""

import numpy as np

# Set up Holland storm field functions
def coriolis(y):
    # Coriolis parameter, input y in radians
    Omega = 2. * np.pi / (60*60*24.) # Angular speed of Earth
    return 2. * Omega * np.sin(np.deg2rad(y))

def dist_to_eye(x,y,eye):
    # xy = (x,y)
    # eye = (x_eye,y_eye)
    xy = np.array([x,y])
    return np.linalg.norm(xy-eye)

def direction_to_eye(x,y,eye):
    xy = np.array([x,y])
    dir_to_eye = (eye-xy) / np.linalg.norm(eye-xy)
    return dir_to_eye

def Wxy(x, y, eye, rw, B, Wmax):
    # Wind speed with x and y coords as input
    # rw - radius of max winds
    # Wmax - maximum wind speed
    # B - model parameter
    r = dist_to_eye(x,y,eye)
    rwrB = (rw/r)**B
    f = coriolis(y)
    return np.sqrt(rwrB * Wmax**2 * np.exp(1-rwrB) + (0.5*r*f)**2) - 0.5*r*f

def pressure(r, rw, B, Pc, Pn):
    # Pressure
    # r - distance from eye
    # rw - radius of max winds
    # Pc - central pressure
    # Pn - atmos pressure at infinity
    # B - model parameter
    return Pc + (Pn - Pc)*np.exp(-(rw/r)**B)

# Define geometry
# Note: all calculations assume we're on a plane
lonmin = -10
lonmax = 10
latmin = 10
latmax = 20
dlon = 0.1
dlat = 0.1
deg2m = 1.1e5 # 1 degree is roughly 1.1*10^5 meters

# Define storm parameters
eye = np.array([-4,15])
B = 1.5
rw = 1e4 / deg2m # Radius of maximum winds (deg)
Wmax = 200 / deg2m # maximum wind speed (deg/s)
Pn = 1005. # Nominal pressure (mb)
Pc = 900. # Central pressure (mb)
storm_u = 0. # storm translation speed (deg/s)
storm_v = 0. # storm translation speed (deg/s)

# Make xy grid
lons = np.arange(lonmin,lonmax,dlon)
lats = np.arange(latmin,latmax,dlat)
nx = len(lons)
ny = len(lats)
dists = np.zeros((ny,nx))
P = np.zeros((ny,nx))
speeds = np.zeros((ny,nx))
U = np.zeros((ny,nx))
V = np.zeros((ny,nx))
for i in range(ny):
    y = lats[i]
    for j in range(nx):
        x = lons[j]
        dists[i,j] = dist_to_eye(x, y, eye)
        P[i,j] = pressure(dists[i,j], rw, B, Pc, Pn)
        speeds[i,j] = Wxy(x, y, eye, rw, B, Wmax)
        uhat,vhat = direction_to_eye(x, y, eye)
        U[i,j] = (vhat * speeds[i,j] + abs(speeds[i,j]/Wmax)*storm_u) * deg2m
        V[i,j] = (-uhat * speeds[i,j] + abs(speeds[i,j]/Wmax)*storm_v) * deg2m


# Write Lat & Lon files
mlon, mlat = np.meshgrid(lons,lats)
np.savetxt('test-lon.dat', mlon, fmt='%.4f')
np.savetxt('test-lat.dat', mlat, fmt='%.4f')

# Write files for U, V, and P
num_times = 8
timestamps = 2016062000 + np.arange(num_times)
output_times = np.repeat(timestamps,ny)
U_out = np.column_stack([output_times,np.row_stack([U]*num_times)])
np.savetxt('test-u10.dat', U_out, fmt='%d'+' %.4f'*nx)
V_out = np.column_stack([output_times,np.row_stack([V]*num_times)])
np.savetxt('test-v10.dat', V_out, fmt='%d'+' %.4f'*nx)
P_out = np.column_stack([output_times,np.row_stack([P]*num_times)])
np.savetxt('test-pmsl.dat', P_out, fmt='%d'+' %.4f'*nx)

# Write trivial bathymetry file
no_data_value = -99999
Z_format = "%7i"
with open("test-bathy.asc", 'w') as outfile:
    # Write out header
    outfile.write('%6i                              ncols\n' % nx)
    outfile.write('%6i                              nrows\n' % ny)
    outfile.write('%22.15e              xlower\n' % lonmin)
    outfile.write('%22.15e              ylower\n' % latmin)
    # write both dx and dy:
    outfile.write('%22.15e    %22.15e          cellsize\n' % (dlon, dlat))
    outfile.write('%10i                          nodata_value\n' % no_data_value)

    # Write out topography data
    Z_format = Z_format + " "
    for i in xrange(ny):
        for j in xrange(nx):
            outfile.write(Z_format % (-5))
        outfile.write("\n")


