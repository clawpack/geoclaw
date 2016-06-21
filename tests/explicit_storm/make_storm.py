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

# Define geometry
# Note: all calculations assume we're on a plane
lonmin = -10
lonmax = 10
latmin = 10
latmax = 20
dlon = 0.05
dlat = 0.05
deg2m = 1.1e5 # 1 degree is roughly 1.1*10^5 meters

# Start time (format is YEARMODAHR)
t0 = 2016062000
# Number of storm snapshots (one per hour)
num_hours = 8

# Define storm parameters
eye = np.array([-4,15])
B = 1.5
rw = 1e4 # Radius of maximum winds (m)
Wmax = 200 # maximum wind speed (m/s)
Pn = 1005. # Nominal pressure (mb)
Pc = 900. # Central pressure (mb)
storm_u = 0. # storm translation speed (m/s)
storm_v = 0. # storm translation speed (m/s)

# Set up Holland storm field functions
def coriolis(y):
    # Coriolis parameter, input y in degrees
    Omega = 2. * np.pi / (60*60*24.) # Angular speed of Earth
    return 2. * Omega * np.sin(np.deg2rad(y))

def dist_to_eye(x,y,eye):
    # eye = (x_eye,y_eye)
    xy = np.array([x,y])
    return np.linalg.norm(xy-eye) * deg2m

def direction_to_eye(x,y,eye):
    xy = np.array([x,y])
    dir_to_eye = (eye-xy) / np.linalg.norm(eye-xy)
    return dir_to_eye

def W(x, y, eye, rw, B, Wmax):
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

# Make xy grid
lons = np.arange(lonmin,lonmax,dlon)
lats = np.arange(latmin,latmax,dlat)
nx = len(lons)
ny = len(lats)
# Create arrays
dists = np.zeros((ny,nx))
P = np.zeros((ny,nx))
speeds = np.zeros((ny,nx))
U = np.zeros((ny,nx))
V = np.zeros((ny,nx))

# Compute storm field
for i in range(ny):
    y = lats[i]
    for j in range(nx):
        x = lons[j]
        dists[i,j] = dist_to_eye(x, y, eye)
        P[i,j] = pressure(dists[i,j], rw, B, Pc, Pn)
        speeds[i,j] = W(x, y, eye, rw, B, Wmax)
        uhat,vhat = direction_to_eye(x, y, eye)
        U[i,j] = vhat * speeds[i,j] + abs(speeds[i,j]/Wmax)*storm_u
        V[i,j] = -uhat * speeds[i,j] + abs(speeds[i,j]/Wmax)*storm_v

# Write Lat & Lon files
mlon, mlat = np.meshgrid(lons,lats)
np.savetxt('test-lon.dat', mlon, fmt='%.4f')
np.savetxt('test-lat.dat', mlat, fmt='%.4f')

# Write files for U, V, and P
timestamps = t0 + np.arange(num_hours)
output_times = np.repeat(timestamps,ny)
U_out = np.column_stack([output_times,np.row_stack([U]*num_hours)])
np.savetxt('test-u10.dat', U_out, fmt='%d'+' %.4f'*nx)
V_out = np.column_stack([output_times,np.row_stack([V]*num_hours)])
np.savetxt('test-v10.dat', V_out, fmt='%d'+' %.4f'*nx)
P_out = np.column_stack([output_times,np.row_stack([P]*num_hours)])
np.savetxt('test-pmsl.dat', P_out, fmt='%d'+' %.4f'*nx)

# Write a trivial bathymetry file
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


