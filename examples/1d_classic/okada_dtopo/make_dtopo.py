"""
Create a dtopo file for 1d geoclaw.
In this example a single fault is split into 3 subfaults for illustration.
If fault.rupture_type == 'kinematic' is used then all subfaults have the
same rupture time and rise time in this example.

The final seafloor deformation with fault.rupture_type == 'kinematic'
should be the same as if fault.rupture_type == 'static' is set, in which
case all subfaults rupture simultaneously and instantaneously at 1 second.
"""

from pylab import *

from clawpack.geoclaw.data import LAT2METER
from importlib import reload

from clawpack.geoclaw import dtopotools
reload(dtopotools)

fault = dtopotools.Fault1d()
fault.subfaults = []

strike = 0.  # 0 ==> top is at right in 1d,  180 ==> top at left
             # assuming 0 < dip <= 90

fault_top_meters = -100e3  # location x0 (in meters) at top edge in 1d
longitude0 = fault_top_meters/LAT2METER  # convert to degrees at latitude 0

fault_top_depth = 20e3     # depth of top below surface (meters)
width = 50e3  # width of fault in 1d (down-dip direction)

dip = 10.  # dip in degrees
theta = dip*pi/180.   # dip in radians


if 0:
    fault.rupture_type = 'static'  # instantaneous and simultaneous on subfaults
    rupture_time = 0.0
    rise_time = 0.
else:
    fault.rupture_type = 'kinematic'
    # use these same values for all subfaults in this simple example:
    rupture_time = 0.0
    rise_time = 20.  
      
average_slip = 10.0

# split into subfaults if desired:
nsubfaults = 3
max_slip = 2*average_slip # if modulated by cosine hump below
dlongitude = width*cos(theta)/LAT2METER / nsubfaults

ddepth = width*sin(theta) / nsubfaults
subfault_width = width/nsubfaults

total_slip = 0.  # keep track
for i in range(nsubfaults):
    # split total slip between subfaults, starting at top
    subfault = dtopotools.SubFault1d()
    subfault.dip = dip
    subfault.width = subfault_width
    subfault.depth = fault_top_depth + ddepth*i
    #subfault.slip = max_slip * 0.5*(1 - cos(2*pi*(i+0.5)/nsubfaults))
    subfault.slip = average_slip  # for constant slip
    total_slip += subfault.slip
    print('subfault %2i at depth %8.3f km has slip = %6.3f' \
            % (i,subfault.depth/1e3,subfault.slip))

    subfault.longitude = longitude0 + i*dlongitude
    subfault.coordinate_specification = 'top center'
    subfault.strike = strike
    subfault.rupture_time = rupture_time # all at same time
    subfault.rise_time = rise_time
    fault.subfaults.append(subfault)

print('average slip = %6.3f' % (total_slip/nsubfaults))


if fault.rupture_type == 'static':
    tend = 1.
    times = [0.,tend]
elif fault.rupture_type == 'kinematic':
    tend = 0.
    for s in fault.subfaults:
        tend = max(tend, s.rupture_time + 2*s.rise_time)
    times = linspace(0,tend,11)
    

print('dtopofile will have times: ',times)
xgrid,zgrid = loadtxt('celledges.data', skiprows=1, unpack=True)

# coarsen if desired:
xgrid = linspace(xgrid[0],xgrid[-1],100)

x = xgrid / LAT2METER  # convert meters to longitude
y = array([0,1])  # for 1d Okada

dtopo2d = fault.create_dtopography(x,y,times)
dtopo = dtopotools.DTopography1d()
dtopo.x = dtopo2d.x * LAT2METER # convert x back from degrees to meters
dtopo.times = dtopo2d.times
dtopo.dZ = dtopo2d.dZ[:,0,:]  # should be constant in y, remove that index

fname = 'dtopo_okada.dtt1'
dtopo.write(fname, dtopo_type=1)
print('Created ',fname)


if 1:
    # make plots
    figure(351,figsize=(6,8))
    clf()
    #dz = dtopo.dZ[-1,:]  # slice in x at final time
    #plot(dtopo.x,dz)
    subplot(211)
    for k in range(1,len(times)):
        plot(dtopo.x/1e3, dtopo.dZ[k,:],label='t = %6.2fs' % times[k])
    title('Okada seafloor deformation')
    legend(loc='upper right',fontsize=6)
    xticks(rotation=20)
    xlim(dtopo.x[0]/1e3,dtopo.x[-1]/1e3)
    xlabel('kilometers')
    ylabel('meters')
    grid(linewidth=0.5)
    
    subplot(212)
    for k,subfault in enumerate(fault.subfaults):
        xtop = subfault.longitude * LAT2METER
        xbottom = xtop + subfault.width*cos(subfault.dip*pi/180)
        ytop = -subfault.depth
        ybottom = ytop - subfault.width*sin(subfault.dip*pi/180)
        plot([xtop/1e3,xbottom/1e3],[ytop/1e3,ybottom/1e3], 
             label='subfault %i' % k)
    title('Subfaults')
    legend(loc='upper right')
    xlim(dtopo.x[0]/1e3,dtopo.x[-1]/1e3)
    xticks(rotation=20)
    xlabel('kilometers')
    ylabel('kilometers below seafloor')
    ylim(1.1*ybottom/1e3,0)
    grid(linewidth=0.5)
    tight_layout()

    fname = 'dtopo_okada.png'
    savefig(fname)
    print('Created ',fname)
    
