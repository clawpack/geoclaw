#!/usr/bin/env python
# coding: utf-8

# # Make Input Files 
# 
# ### For `$CLAW/geoclaw/examples/tsunami/eta_init_force_dry`
# 
# For this example simple artificial topography is generated in order to illustrate various things.
# 
# Contents:
# 
#  - [Define ocean topography](#topo_ocean)
#  - [Define topo for small coastal region](#topo_coast)
#  - [Create dtopo for an earthquake source](#dtopo)
#  - [Force Dry array](#force_dry)
#  
# Running this notebook should create a set of files in the directory `input_files`.
# 
# Alternatively, running 
# 
#     make input
#     
# or equivalently 
# 
#     python make_input_files.py
#     
# will run the python script version of this notebook, which was created with the command
# 
#     jupyter nbconvert --to python --TagRemovePreprocessor.enabled=True \
#             --TagRemovePreprocessor.remove_cell_tags="['hide-py']" \
#             make_input_files.ipynb
# 
# This will only work if [nbconvert](https://nbconvert.readthedocs.io/en/latest/index.html) is installed.
# 
# Note that cells in this notebook that create plots are not included in the `.py` version (due to the cell tag `hide-py` that is applied to these cells, visible if you select `View -> Cell Toolbar -> Tags` in the notebook menu).

# In[2]:


from pylab import *
from scipy.interpolate import interp1d
import os


# In[3]:


from clawpack.geoclaw import topotools, marching_front, dtopotools
from clawpack.visclaw import plottools


# ## Directory for input files:

# In[4]:


inputdir = 'input_files'
os.system('mkdir -p %s' % inputdir)
print('Input files will be put in directory %s' % inputdir)


# <div id="topo_ocean"></div>
# 
# ## Define ocean topography
# 
# This simple topography is piecewise linear in $x$ (longitude) with a continental shelf and beach, and constant in the $y$ (latitude) direction.  It is placed at the equator so distances are roughly equal in $x$ and $y$, and also placed at longitude 0.  

# In[5]:


# Define piecewise linear function (unequally spaced):
xocean = array([-2,-1,-0.5,-0.1,0.1])
zocean = array([-3000,-3000,-100,-100,100])

# Interpolate to equally spaced grid for topofile:
xo = arange(-2,0.2,0.1)
yo = array([-2,2])
zfunc = interp1d(xocean,zocean,fill_value="extrapolate")
zo = zfunc(xo)

# Convert to 2d arrays:
Xo,Yo = meshgrid(xo,yo)
Zo = vstack((zo,zo))


# ### Save as a topofile:

# In[7]:


topo = topotools.Topography()
topo.set_xyZ(xo,yo,Zo)

topofile = '%s/topo_ocean.tt3' % inputdir
topo.write(topofile, topo_type=3, Z_format="%11.3e")
print('Created ', topofile)


# <div id="topo_coast"></div>
# 
# ## Define topo for small coastal region
# 
# We define some more complicated topography on a finer grid over a small coastal region with 1/3 arcsecond resolution, chosen to be aligned with integer multiples of degrees (e.g. a grid point at longitude `x=0` and latitude `y=0`) as typical of real DEMs from NCEI.  This is important when aligning computational grids and fgmax grids (if used) in `setrun.py`.   
# 
# We will use a cutoff function so that this fine-scale topo matches the linear beach profile of the ocean topography along the edges of this rectangle.  The cutoff is 1 in the center of the rectangle and decays to 0 at the edges:

# In[8]:


# choose DEM grid points:
arcsec13 = 1./(3*3600.)  # 1/3 arcsecond
print('arcsec13 = %.6f degrees = %.2f meters' % (arcsec13,arcsec13*111e3))
x = arange(-100*arcsec13, 150*arcsec13, arcsec13)
y = arange(-55*arcsec13, 55*arcsec13, arcsec13)
X,Y = meshgrid(x,y)
print('X.shape = ', X.shape)

x1,x2 = x.min(), x.max()
y1,y2 = y.min(), y.max()
print('Extent of coastal topo: (%.6f, %.6f, %.6f, %.6f)' % (x1,x2,y1,y2))

# define the cutoff function:

w = 0.001 # width of cutoff layer
cutoff = 1. / (1. + exp(1e4*(X-(x2-w))) + exp(1e4*((x1+w)-X))                   + exp(1e4*(Y-(y2-w))) + exp(1e4*((y1+w)-Y)))


# The topography in this region is the linearly sloping beach augmented by a Gaussian dip.  The beach slope is chosen to agree with the ocean topography offshore (1 km / degree, about 1/100), while onshore there is a smaller slope in this region for illustration.

# In[10]:


Z0 = 1e3*X  # sloping beach matching ocean topography
Z1 = where(X<0, 1e3*X, 0.2e3*X)  # smaller slope on shore
R1 = (X-0.004)**2 + (Y-0.002)**2
Z1 += -4*exp(-500000*R1)         # Gaussian dip
Z = (1-cutoff)*Z0 + cutoff*Z1


# The lower plot in the figure above shows the same topography as on the top, but with x,y units of meters to better show the scale.  Recall that 1 degree is about 111 km and 1/3 arcsec is about 10 meters.
# 
# In the plots above, the red contour is at $Z = 0$, and hence is the "shoreline".  However, the isolated "lake" with elevation $Z < 0$ could be dry land below sea level.  Normally with GeoClaw this region would be filled with water initially up to $Z = 0$ everywhere.  Below in [the Force_Dry section](#force_dry), we discuss how to force this region to be initialized as dry if it is in fact dry land.

# ### Save this as a topofile:

# In[12]:


topo = topotools.Topography()
topo.set_xyZ(x,y,Z)

topofile = '%s/topo_shore.tt3' % inputdir
topo.write(topofile, topo_type=3, Z_format="%11.3e")
print('Created ', topofile)


# In the plot on the left above, the black rectangle showing the extent of the coastal DEM is barely visible.  Zooming in shows that the topography does match up near the edges of this rectangle.  In GeoClaw the finest available topography is used when computing cell-averaged topo values, so the coastal DEM will be used for any cell that overlaps this region. 

# <div id="dtopo"></div>
# 
# ## Create dtopo for an earthquake source:
# 
# We define a simple earthquake in which there is uniform slip on a single subfault. The parameters are chosen to be somewhat reasonable for a subduction zone event offshore, but the shape is a bit odd (width 100 km and length 50 km) in order to give a smallish event with the desired onshore subsidence, for illustration purposes.

# In[14]:


subfault = dtopotools.SubFault()
subfault.strike = 0.
subfault.length = 50.e3
subfault.width = 100.e3
subfault.depth = 10.e3
subfault.slip = 5.
subfault.rake = 90.
subfault.dip = 10.
subfault.longitude = -1.
subfault.latitude = 0.
subfault.coordinate_specification = "top center"

fault = dtopotools.Fault()
fault.subfaults = [subfault]

print("Earthquake magnitude: Mw = %.2f" % fault.Mw())
dtopo_fname = '%s/dtopo_test.tt3' % inputdir
print("Using Okada model to create dtopo file", dtopo_fname)

x_deform = linspace(-2, 1, 100)
y_deform = linspace(-1, 1, 100)
times = [1.]

fault.create_dtopography(x_deform,y_deform,times)
dtopo = fault.dtopo

dtopo.write(dtopo_fname, dtopo_type=3)


# The left plot above shows the sea floor deformation as contours and colors, along with the extent of the continental shelf as blue dashed lines and the shoreline as a red dashed line. The plot on the right shows the vertical deformation along a transect at latitude 0 going through the coastal region of interest.  
# 
# We can compute the subsidence at the location on the shoreline where our fine scale topography is defined as:

# In[16]:


xlon = 0.
ilon = where(dtopo.x<=xlon)[0].max()
ylat = 0.
jlat = where(dtopo.y<=ylat)[0].max()
#print(ilon,jlat)
dz0 = dtopo.dZ[0,jlat,ilon]
print('Surface deformation at x=%.2f, y=%.2f is dz = %.2f meters'       % (xlon,ylat,dz0))


# <div id="force_dry"></div>
# 
# # Force Dry array
# 
# Now suppose that the onshore lake shown in the plots above is really a depression that should be dry land in spite of being below sea level.  We can use the marching front algorithm from [`clawpack.geoclaw.marching_front`](http://depts.washington.edu/clawpack/sampledocs/dev_v5.7.0/marching_front.html) to identify points that are below sea level but disconnected from the coast.  
# 
# We use the marching front algorithm starting by assuming any point with `Z < Z1 = -5` meters should be wet and marching to find all connected points with elevation up to `Z = Z2 = 0`:

# In[18]:


wet_points = marching_front.select_by_flooding(topo.Z, Z1=-5., Z2=0., max_iters=None)


# See the documentation page [Force Cells to be Dry Initially](http://depts.washington.edu/clawpack/sampledocs/dev_v5.7.0/force_dry.html) for more discussion of the cells below...

# ## Create `force_dry_init` array for GeoClaw
# 
# First we buffer the points identified above as discussed in the the documentation page [Force Cells to be Dry Initially](http://depts.washington.edu/clawpack/sampledocs/dev_v5.7.0/force_dry.html).

# In[20]:


dry_points = 1 - wet_points
dry_points_sum = dry_points[1:-1,1:-1] + dry_points[0:-2,1:-1] + dry_points[2:,1:-1] +                  dry_points[1:-1,0:-2] + dry_points[0:-2,0:-2] + dry_points[2:,0:-2] +                  dry_points[1:-1,2:] + dry_points[0:-2,2:] + dry_points[2:,2:]
        
# initialize array to 0 everywhere:
force_dry_init = zeros(dry_points.shape)

# reset in interior to 1 if all points in the 3x3 block around it are dry:
force_dry_init[1:-1,1:-1] = where(dry_points_sum == 9, 1, 0)


# And finally create the input file needed for GeoClaw.  Note that this creates a file with the same format as a topofile having `topo_type == 3` as described in [Topography Data documentation](http://www.clawpack.org/topo.html).  We specify `Z_format= '%1i'` to print out single-digit integers since this file has values 0 or 1 rather than topography elevations (with 1 indicated points that should be forced to be dry when initializing grid patches in GeoClaw).

# In[21]:


force_dry_init_topo = topotools.Topography()
force_dry_init_topo.set_xyZ(topo.x,topo.y,force_dry_init)

fname_force_dry_init = '%s/force_dry_init.tt3' % inputdir
force_dry_init_topo.write(fname_force_dry_init, topo_type=3, Z_format='%1i')
print('Created %s' % fname_force_dry_init)


# See [run_geoclaw.ipynb](run_geoclaw.ipynb) for more discussion and sample GeoClaw results.

# In[ ]:




