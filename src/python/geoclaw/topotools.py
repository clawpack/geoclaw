#!/usr/bin/env python
# encoding: utf-8

r"""GeoClaw Topography Tools Module

Module provides several functions for reading, writing and manipulating
topography (bathymetry) files.

:Functions:
 - dms2decimal - Convert (degrees, minutes, seconds) to decimal degrees
 - dist_meters2latlong - Convert dx, dy distance in meters to degrees
 - dist_latlong2meters - Convert dx, dy distance in degrees to meters
 - haversine - Calculate the haversine based great circle distance
 - inv_haversine - Inverts the haversine distance

:Topography Class:

:TODO:
 - Tests are implemented but not passing, should we expect the arrays to be 
   identical?
 - Add sub and super sampling capababilities
 - Add functions for creating topography based off a topo function, incorporate
   the create_topo_func into Topography class, maybe allow more broad 
   initialization ability to the class to handle this?
 - Fix `in_poly` function
 - Add remove/fill no data value
 - Probably should better handle remote files (fetching from http)
 - Add more robust plotting capabilities
"""

import os
import urllib
import types

import numpy



# Constants
from data import Rearth
DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi


# ==============================================================================
#  Functions
# ==============================================================================
def dms2decimal(d,m,s,coord='N'):
    r"""Convert coordinates in (degrees, minutes, seconds) to decimal form.  
    
    If coord == 'S' or coord == 'W' then value is negated too.
    Example: 
        >>> topotools.dms2decimal(7,30,36,'W')
        -7.51
    (Note that you might want to add 360 to resulting W coordinate
    if using E coordinates everywhere in a computation spanning date line.)

    returns float

    """

    deg = d + m / 60.0 + s / 3600.0
    if coord in ['S','W']:
        deg = -deg

    return deg


def dist_latlong2meters(dx, dy, latitude=0.0):
    """Convert distance from degrees longitude-latitude to meters.

    Takes the distance described by *dx* and *dy* in degrees and converts it into
    distances in meters.

    returns (float, float) 

    """

    dym = Rearth * DEG2RAD * dy
    dxm = Rearth * numpy.cos(latitude * DEG2RAD) * dx * DEG2RAD

    return dxm,dym


def dist_meters2latlong(dx, dy, latitude=0.0):
    """Convert distance from meters to degrees of longitude-latitude.

    Takes the distance described by *dx* and *dy* in meters and converts it into
    distances in the longitudinal and latitudinal directions in degrees.  

    returns (float, float)

    """

    dxd = dx / (Rearth * numpy.cos(latitude * DEG2RAD)) * RAD2DEG
    dyd = dy * RAD2DEG / Rearth

    return dxd, dyd


def haversine(x, y, units='degrees'):
    r"""Compute the great circle distance on the earth between points x and y.


    """
    if units == 'degrees':
        # convert to radians:
        x *= DEG2RAD
        y *= DEG2RAD

    delta = [x[0] - y[0], x[1] - y[1]]

    # angle subtended by two points, using Haversine formula:
    dsigma = 2.0 * numpy.arcsin( numpy.sqrt( numpy.sin(0.5 * delta[1])**2   \
            + numpy.cos(x[1]) * numpy.cos(y[1]) * numpy.sin(0.5 * delta[0])**2))

    # alternative formula that may have more rounding error:
    #dsigma2 = arccos(sin(y1)*sin(y2)+ cos(y1)*cos(y2)*cos(dx))
    #print "max diff in dsigma: ", abs(dsigma-dsigma2).max()

    return Rearth * dsigma
    

def inv_haversine(d,x1,y1,y2,Rsphere=Rearth,units='degrees'):
    r"""Invert the Haversine function to find dx given a distance and point.


    Invert the haversine function to find dx given distance d and (x1,y1) and y2.
    The corresponding x2 can be x1+dx or x1-dx.
    May return NaN if no solution.
    """

    if units=='degrees':
        # convert to radians:
        x1 = x1 * RAD2DEG
        y1 = y1 * RAD2DEG
        y2 = y2 * RAD2DEG
    elif units != 'radians':
        raise Exception("unrecognized units")
    dsigma = d / Rsphere
    cos_dsigma = (numpy.cos(dsigma) - numpy.sin(y1)*numpy.sin(y2)) / (numpy.cos(y1)*numpy.cos(y2))
    dx = numpy.arccos(cos_dsigma)
    if units=='degrees':
        dx = dx * RAD2DEG
    return dx


def create_topo_func(loc,verbose=False):
    """
    Given a 1-dimensional topography profile specfied by a set of (x,z) 
    values, create a lambda function that when evaluated will give the 
    topgraphy at the point (x,y).  (The resulting function is constant in y.)
    
    :Example: 
    >>> f = create_topo_profile_func(loc)
    >>> b = f(x,y)
    
    :Input:
     - *loc* (list) - Create a topography file with the profile denoted by the
       tuples inside of loc.  A sample set of points are shown below.  Note 
       that the first value of the list is the x location and the second is 
       the height of the topography.
        
        z (m)
        ^                                                  o loc[5]  o
        |                                                    
        |                                          loc[4]   
        |--------------------------------------------o-----> x (m) (sea level)
        |                                            
        |                                o loc[2] o loc[3]
        |                         
        |                         
        |                           o loc[1]
        |           
        |                               
        |__________________o loc[0]
        0.0               
        
        
    """
    
    cmd_str = "lambda x,y: (x <= %s) * %s" % (loc[0][0],loc[0][1])
    for i in xrange(0,len(loc)-1):
        loc_str = " + (%s < x) * (x <= %s)" % (loc[i][0],loc[i+1][0])
        loc_str = "".join((loc_str," * ((%s - %s) " % (loc[i][1],loc[i+1][1])))
        loc_str = "".join((loc_str," / (%s - %s)" % (loc[i][0],loc[i+1][0])))
        loc_str = "".join((loc_str," * (x - %s) + %s)" % (loc[i][0],loc[i][1])))
        cmd_str = "".join((cmd_str,loc_str))
    cmd_str = "".join((cmd_str," + (%s < x) * %s" % (loc[-1][0],loc[-1][1])))
    
    if verbose:
        print cmd_str
    return eval(cmd_str)


def topo1writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints):
    """
    Function topo1writer will write out the topofiles by evaluating the
    function topo on the grid specified by the other parameters.

    Assumes topo can be called on arrays X,Y produced by numpy.meshgrid.

    Output file is of "topotype1," which we use to refer to a file with
    (x,y,z) values on each line, progressing from upper left corner across
    rows, then down.
    """
    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)
    
    topography.write(outfile, topo_type=1)


def topo2writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
    r"""Write out a topo type 2 file by evaluating the function *topo*.

    This routine is here for backwards compatibility and simply creates a new
    topography object and writes it out.

    """

    topography = Topography(topo_func=topo)

    topography.x = numpy.linspace(xlower,xupper,nxpoints)
    topography.y = numpy.linspace(ylower,yupper,nypoints)

    topography.write(outfile, topo_type=2)


def get_topo(topo_fname, remote_directory, force=None):
    """
    Download a topo file from the web, provided the file does not
    already exist locally.

    remote_directory should be a URL.  For GeoClaw data it may be a
    subdirectory of  http://kingkong.amath.washington.edu/topo/
    See that website for a list of archived topo datasets.

    If force==False then prompt the user to make sure it's ok to download,
    with option to first get small file of metadata.

    If force==None then check for environment variable CLAW_TOPO_DOWNLOAD
    and if this exists use its value.  This is useful for the script
    python/run_examples.py that runs all examples so it won't stop to prompt.
    """
    import urllib

    if force is None:
        CTD = os.environ.get('CLAW_TOPO_DOWNLOAD', None)
        force = (CTD in [True, 'True'])
    print 'force = ',force

    if os.path.exists(topo_fname):
        print "*** Not downloading topo file (already exists): %s " % topo_fname
    else:
        remote_fname = topo_fname
        local_fname = topo_fname
        remote_fname_txt = remote_fname + '.txt'
        local_fname_txt = local_fname + '.txt'

        print "Require remote file ", remote_fname
        print "      from ", remote_directory
        if not force:
            ans=raw_input("  Ok to download topo file?  \n"  +\
                          "     Type y[es], n[o] or ? to first retrieve and print metadata  ")
            if ans.lower() not in ['y','yes','?']:
                print "*** Aborting!   Missing: ", local_fname
                return
            if ans=="?":
                try:
                    print "Retrieving remote file ", remote_fname_txt
                    print "      from ", remote_directory
                    url = os.path.join(remote_directory, remote_fname_txt)
                    urllib.urlretrieve(url, local_fname_txt)
                    os.system("cat %s" % local_fname_txt)
                except:
                    print "*** Error retrieving metadata file!"
                ans=raw_input("  Ok to download topo file?  ")
                if ans.lower() not in ['y','yes','?']:
                    print "*** Aborting!   Missing: ", local_fname
                    return

        if not os.path.exists(local_fname_txt):
            try:
                print "Retrieving metadata file ", remote_fname_txt
                print "      from ", remote_directory
                url = os.path.join(remote_directory, remote_fname_txt)
                urllib.urlretrieve(url, local_fname_txt)
            except:
                print "*** Error retrieving metadata file!"

        try:
            print "Retrieving topo file ", remote_fname
            print "      from ", remote_directory
            url = os.path.join(remote_directory, remote_fname)
            urllib.urlretrieve(url, local_fname)
        except:
            print "*** Error retrieving file!  Missing: ", local_fname
            raise Exception("Error from urllib.urlretrieve")
        try:
            firstline = open(local_fname,'r').readline()
            if firstline.find('DOC') > -1:
                print "*** Possible error -- check the file ", local_fname
            else:
                print "Saved to ", local_fname
        except:
            raise Exception("Error opening file %s" % local_fname)


def swapheader(inputfile, outputfile):
    r"""Swap the order of key and value in header to value first.

    Note that this is a wrapper around functionality in the Topography class.

    """
    topo = Topography(inputfile)
    topo.write(outputfile)


# ==============================================================================
#  Topography class
# ==============================================================================
class Topography(object):

    r"""Base topography class.

    A class representing a single topography file.

    :Properties:

    :Initialization:
         - 

    :Examples:

        >>> import clawpack.geoclaw.topotools as topo
        >>> topo_file = topo.Topography('./topo.tt3')
        >>> topo_file.plot()

    """

    @property
    def z(self):
        r"""A representation of the data as an 1d array."""
        if self._z is None:
            self.read(mask=True)
        return self._z
    @z.setter
    def z(self, value):
        self._z = value
    @z.deleter
    def z(self):
        del self._z

    @property
    def Z(self):
        r"""A representation of the data as a 2d array."""
        if self._Z is None:
            self.generate_2d_depths(mask=True)
        return self._Z
    @Z.setter
    def Z(self, value):
        self._Z = value
    @Z.deleter
    def Z(self):
        del self._Z

    @property
    def x(self):
        r"""One dimensional coorindate array in x direction."""
        if self._x is None:
            self.read(mask=True)
        return self._x
    @x.setter
    def x(self, value):
        self._extent = None
        self._x = value
    @x.deleter
    def x(self):
        del self._x

    @property
    def X(self):
        r"""Two dimensional coordinate array in x direction."""
        if self._X is None:
            self.generate_2d_coordinates(mask=True)
        return self._X
    @X.deleter
    def X(self):
        del self._X

    @property
    def y(self):
        r"""One dimensional coordinate array in y direction."""
        if self._y is None:
            self.read(mask=True)
        return self._y
    @y.setter
    def y(self, value):
        self._extent = None
        self._y = value
    @y.deleter
    def y(self):
        del self._y

    @property
    def Y(self):
        r"""Two dimensional coordinate array in y direction."""
        if self._Y is None:
            self.generate_2d_coordinates(mask=True)
        return self._Y
    @Y.deleter
    def Y(self):
        del self._Y

    @property
    def extent(self):
        r"""Extent of the topography."""
        if self._extent is None:
            self._extent = ( numpy.min(self.x), numpy.max(self.x), 
                             numpy.min(self.y), numpy.max(self.y) )
        return self._extent
    @extent.setter
    def extent(self, value):
        self._extent = value

    @property
    def delta(self):
        r"""Spacing of data points."""
        if self._delta is None:
            if self.unstructured:

                # Calculate the smallest spacing between grid points            
                dx = numpy.infty
                dy = numpy.infty
                num_comparisons = self.x.shape[0] - 1
                for i in xrange(self.x.shape[0]):
                    for j in xrange(num_comparisons):
                        dx = min(dx, numpy.abs(self.x[i + j + 1] - self.x[i]))
                        dy = min(dy, numpy.abs(self.y[i + j + 1] - self.y[i]))

                    num_comparisons -= 1
                self._delta = [dx, dy]
            else:
                # All other topography types should have equally spaced grid
                # points in each direction
                begin_delta = numpy.array([abs(self.x[1] - self.x[0]), 
                                           abs(self.y[1] - self.y[0])])
                end_delta =   numpy.array([abs(self.x[-2] - self.x[-1]), 
                                           abs(self.y[-2] - self.y[-1])])
                assert numpy.allclose(begin_delta, end_delta, 1e-8),   \
                       "Grid spacing delta not constant, %s != %s." %  \
                       (begin_delta, end_delta)
                self._delta = numpy.round(begin_delta[0], 15) 
        return self._delta


    def __init__(self, path=None, topo_func=None, topo_type=None, 
                       unstructured=False, force=False):
        r"""Topography initialization routine.
        
        See :class:`Topography` for more info.

        """

        super(Topography, self).__init__()

        self.path = path
        self.topo_func = topo_func
        if self.path is None:
            # We are going to generate the topography via the provided
            # topo_func function
            if (topo_func is None or 
                not isinstance(topo_func, types.FunctionType)):
                raise ValueError("Must provide either a path to a topography ",
                                 "file or a generator function.")

            # Do nothing for right now, wait until user fills in data
        else:
            if topo_type is not None:
                self.topo_type = topo_type
            else:
                # Try to look at suffix for type
                extension = os.path.splitext(path)[1][1:]
                if extension[:2] == "tt":
                    self.topo_type = int(extension[2])
                elif extension == 'xyz':
                    self.topo_type = 1
                elif extension == 'asc':
                    self.topo_type = 3
                else:
                    # Default to 3
                    self.topo_type = 3

            # Check if the path is a URL and fetch data if needed or forced
            if "http" in self.path:
                new_path = os.path.join(os.getcwd(), os.path.split(self.path)[0])
                if not os.path.exists(new_path) or force:
                    urllib.urlretrieve(self.path)

                # Change path to be local
                self.path = new_path

        self.unstructured = unstructured
        self.no_data_value = -9999

        # Data storage for only calculating array shapes when needed
        self._z = None
        self._Z = None
        self._x = None
        self._X = None
        self._y = None
        self._Y = None
        self._extent = None
        self._delta = None

        self.coordinate_transform = lambda x,y: (x,y)


    def generate_2d_depths(self, mask=True):
        r"""Generate a 2d array of the depths."""

        # Check to see if we need to generate these
        if self._Z is None:

            if self.unstructured:
                # Really no way to do this here with performing a projection via
                # extract.  Note that if the projection is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d arrays, first project the data and" \
                                 + " try to perform this operation again.") 

            if self.path is not None:
                if self._z is None:
                # Try to read the data, may not have done this yet
                    self.read(mask=mask)
                    if self._Z is not None:
                        # We are done, the read function did our work
                        return

                # See if self._X and self._Y are already computed and use them if
                # available, otherwise just use self._x and self._y
                if self._X is not None and self._Y is not None:
                    new_shape = self._X.shape
                else:
                    new_shape = (self._x.shape[0], self._y.shape[0])
                # Reshape, note that the mask follows along with the new array
                self._Z = numpy.reshape(self._z, new_shape)

            elif self.topo_func is not None:
                # Generate topo via topo_func
                self._Z = numpy.flipud(self.topo_func(self.X, self.Y))


    def generate_2d_coordinates(self, mask=True):
        r"""Generate 2d coordinate arrays."""

        # Check to see if we need to generate these
        if self._X is None and self._Y is None:

            if self.unstructured:
                # Really no way to do this here with performing a projection via
                # extract.  Note that if the projection is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d coordinates, first project the data" \
                                 + " and try to perform this operation again.")

            if self.path is not None:
                if abs(self.topo_type) == 1:
                    # Reading this topo_type should produce the X and Y arrays
                    self.read(mask=mask)
                elif abs(self.topo_type) in [2,3]:
                    if self._x is None or self._y is None:
                        # Try to read the data to get these, may not have been done yet
                        self.read(mask=mask)
                    # Generate arrays
                    self._X, self._Y = numpy.meshgrid(self._x, self._y)
                else:
                    raise ValueError("Unrecognized topo_type: %s" % self.topo_type)

            elif self.topo_func is not None:
                if self._x is None or self._y is None:
                    raise ValueError("The x and y arrays must be set to ",
                                     "create 2d coordinate arrays.")
                self._X, self._Y = numpy.meshgrid(self._x, self._y)

            
            # If masking has been requested try to get the mask first from 
            # self._Z and then self._z
            if mask:
                if self._Z is None:
                    # Check to see if we really need to do anything here
                    if isinstance(self._z, numpy.ma.MaskedArray):
                        # Try to create self._Z
                        self.generate_2d_depths(mask=mask)

                if isinstance(self._Z, numpy.ma.MaskedArray):
                    # Use Z's mask for the X and Y coordinates
                    self._X = numpy.ma.MaskedArray(self._X, mask=self._Z.mask, 
                                                                     copy=False)
                    self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Z.mask, 
                                                                     copy=False)


    def read(self, mask=True, filter_region=None):
        r"""Read in the data from the object's *path* attribute.

        Stores the resulting data in one of the sets of *x*, *y*, and *z* or 
        *X*, *Y*, and *Z*.  

        :Input:
         - *mask* (bool)
         - *filter_region* (tuple)

        """

        if self.unstructured:
            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)
            points = []
            values = []

            # Filter region if requested
            if filter_region is not None:
                for coordinate in data:
                    if filter_region[0] <= coordinate[0] <= filter_region[1]:
                        if filter_region[2] <= coordinate[1] <= filter_region[3]:
                            points.append(coordinate[0:2])
                            values.append(coordinate[2])

                if len(points) == 0:
                    raise Exception("No points were found inside requested " \
                                  + "filter region.")

                # Cast lists as ndarrays
                self._x = numpy.array(points[:,0])
                self._y = numpy.array(points[:,1])
                self._z = numpy.array(values)

            else:
                self._x = data[:,0]
                self._y = data[:,1]
                self._z = data[:,2]

        else:
            # Data is in one of the GeoClaw supported formats
            if abs(self.topo_type) == 1:
                data = numpy.loadtxt(self.path)
                N = [0,0]
                y0 = data[0,1]
                for (n, y) in enumerate(data[1:,1]):
                    if y != y0:
                        N[1] = n + 1
                        break
                N[0] = data.shape[0] / N[1]

                self._X = data[:,0].reshape(N)
                self._x = data[:N[0],0]
                self._Y = data[:,1].reshape(N)
                self._y = data[::N[0],1]
                self._Z = data[:,2].reshape(N)
                self._delta = self.X[0,1] - self.X[0,0]

            elif abs(self.topo_type) in [2,3]:
                # Get header information
                N = self.read_header()
                self._x = numpy.linspace(self.extent[0], self.extent[1], N[0])
                self._y = numpy.linspace(self.extent[2], self.extent[3], N[1])

                if abs(self.topo_type) == 2:
                    # Data is read in as a single column, reshape it
                    self._Z = numpy.loadtxt(self.path, skiprows=6).reshape(N[1],N[0])
                    self._Z = numpy.flipud(self._Z)
                elif abs(self.topo_type) == 3:
                    # Data is read in starting at the top right corner
                    self._Z = numpy.flipud(numpy.loadtxt(self.path, skiprows=6))
        
                if mask:
                    self._Z = numpy.ma.masked_values(self._Z, self.no_data_value, copy=False)
                    
            else:
                raise IOError("Unrecognized topo_type: %s" % self.topo_type)
                
            if self.topo_type < 0:
                # positive Z is depth below sea level, so negate:
                self._Z = -self._Z
                
                
            # Perform region filtering
            if filter_region is not None:
                # Find indices of region
                region_index = [None, None, None, None]
                region_index[0] = (self.x >= filter_region[0]).nonzero()[0][0]
                region_index[1] = (self.x <= filter_region[1]).nonzero()[0][-1]
                region_index[2] = (self.y >= filter_region[2]).nonzero()[0][0]
                region_index[3] = (self.y <= filter_region[3]).nonzero()[0][-1]

                self._x = self._x[region_index[0]:region_index[1]]
                self._y = self._y[region_index[2]:region_index[3]]

                # Force regeneration of 2d coordinate arrays and extent
                if self._X is not None or self._Y is not None:
                    del self._X, self._Y
                    self._X = None
                    self._Y = None
                self._extent = None

                # Modify Z array as well
                self._Z = self._Z[region_index[2]:region_index[3],
                                  region_index[0]:region_index[1]]


    def read_header(self):
        r"""Read in header of topography file at path.

        If a value returns numpy.nan then the value was not retrievable.  Note
        that this routine can read in headers whose values and labels are 
        swapped.

        """

        if abs(self.topo_type) in [2,3]:

            # Default values to track errors
            num_cells = [numpy.nan,numpy.nan]
            self._extent = [numpy.nan,numpy.nan,numpy.nan,numpy.nan]
            self._delta = numpy.nan

            with open(self.path, 'r') as bathy_file:
                # Check to see if we need to flip the header values
                first_line = bathy_file.readline()
                try:
                    num_cells[0] = int(first_line.split()[0])
                except ValueError:
                    # Assume the header is flipped from what we expect
                    num_cells[0] = int(first_line.split()[-1])
                    value_index = -1
                else:
                    value_index = 0

                num_cells[1] = int(bathy_file.readline().split()[value_index])
                self._extent[0] = float(bathy_file.readline().split()[value_index])
                self._extent[2] = float(bathy_file.readline().split()[value_index])
                self._delta = float(bathy_file.readline().split()[value_index])
                self.no_data_value = float(bathy_file.readline().split()[value_index])
                
                self._extent[1] = self._extent[0] + num_cells[0] * self.delta
                self._extent[3] = self._extent[2] + num_cells[1] * self.delta

        else:
            raise IOError("Cannot read header for topo_type %s" % self.topo_type)
            
        return num_cells

    def write(self, path, no_data_value=None, topo_type=None, masked=True):
        r"""Write out a topography file to path of type *topo_type*.

        Writes out a bathymetry file of topo type specified with *topo_type* or
        inferred from the output file's extension, defaulting to 3, to path
        from data in Z.  The rest of the arguments are used to write the header
        data.

        """

        # Determine topo type if not specified
        if topo_type is None:
            if self.topo_type is not None:
                topo_type = self.topo_type
            else:
                # Try to look at suffix for type
                extension = os.path.splitext(path)[1][1:]
                if extension[:2] == "tt" or extension[:2] == 'topotype':
                    topo_type = int(extension[2])
                elif extension == 'xyz':
                    topo_type = 1
                else:
                    # Default to 3
                    topo_type = 3

        # Default to this object's no_data_value if the passed is None, 
        # otherwise the argument will override the object's value or it will 
        # default to -9999 (default for the class)
        if no_data_value is None:
            no_data_value = self.no_data_value

        # Check to see if masks have been applied to bathymetry, if so use them
        # if masked is True
        if isinstance(self.Z, numpy.ma.MaskedArray) and masked:
            pass
        else:
            pass

        with open(path, 'w') as outfile:
            if self.unstructured:
                for (i, depth) in enumerate(self.z):
                    outfile.write("%s %s %s\n" % (self.x[i], self.y[i], depth))

            elif topo_type == 1:
                # longitudes = numpy.linspace(lower[0], lower[0] + delta * Z.shape[0], Z.shape[0])
                # latitudes = numpy.linspace(lower[1], lower[1] + delta * Z.shape[1], Z.shape[1])
                for (j, latitude) in enumerate(self.y):
                    for (i, longitude) in enumerate(self.x):
                        outfile.write("%s %s %s\n" % (longitude, latitude, self.Z[j,i]))

            elif topo_type == 2 or topo_type == 3:
                # Write out header
                outfile.write('%6i                              ncols\n' % self.Z.shape[1])
                outfile.write('%6i                              nrows\n' % self.Z.shape[0])
                outfile.write('%22.15e              xlower\n' % self.extent[0])
                outfile.write('%22.15e              ylower\n' % self.extent[2])
                outfile.write('%22.15e              cellsize\n' % self.delta)
                outfile.write('%10i                          nodata_value\n' % no_data_value)

                masked_Z = isinstance(self.Z, numpy.ma.MaskedArray)

                # Write out bathy data
                if topo_type == 2:
                    if masked_Z:
                        Z_filled = numpy.flipud(self.Z.filled())
                    else:
                        Z_filled = numpy.flipud(self.Z)
                    for i in xrange(self.Z.shape[0]):
                        for j in xrange(self.Z.shape[1]):
                            outfile.write("%22.15e\n" % Z_filled[i,j])
                    if masked_Z:
                        del Z_filled
                elif topo_type == 3:
                    if masked_Z:
                        Z_flipped = numpy.flipud(self.Z.filled())
                    else:
                        Z_flipped = numpy.flipud(self.Z)
                    for i in xrange(self.Z.shape[0]):
                        for j in xrange(self.Z.shape[1]):
                            outfile.write("%22.15e   " % (Z_flipped[i,j]))
                        outfile.write("\n")
                    if masked_Z:
                        del Z_flipped

            else:
                raise NotImplemented("Output type %s not implemented." % topo_type)


    def plot(self, axes=None, region_extent=None, contours=None, 
             coastlines=True, limits=None, cmap=None, fig_kwargs={}):
        r"""Plot the topography."""

        import matplotlib.pyplot as plt
        import clawpack.visclaw.colormaps as colormaps

        # Create axes if needed
        if axes is None:
            fig = plt.figure(**fig_kwargs)
            axes = fig.add_subplot(111)
        
        # Turn off annoying offset
        axes.ticklabel_format(format="plain", useOffset=False)

        # Generate limits if need be
        if region_extent is None:
            region_extent = ( numpy.min(self.X), numpy.max(self.X),
                              numpy.min(self.Y), numpy.max(self.Y) )
        mean_lat = 0.5 * (region_extent[3] - region_extent[2])
        axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))
        if limits is None:
            depth_extent = (numpy.min(self.Z),numpy.max(self.Z))
        else:
            depth_extent = limits

        # Create color map
        if cmap is None:
            land_cmap = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                                 0.25:[0.0,1.0,0.0],
                                                  0.5:[0.8,1.0,0.5],
                                                  1.0:[0.8,0.5,0.2]})
            sea_cmap = plt.get_cmap('Blues_r')
            cmap = colormaps.add_colormaps((land_cmap, sea_cmap), 
                                           data_limits=depth_extent,
                                           data_break=0.0)

        # Plot data
        if contours is not None:
            plot = axes.contourf(self.X, self.Y, self.Z, contours,cmap=cmap)
        elif isinstance(self.Z, numpy.ma.MaskedArray):
            plot = axes.pcolor(self.X, self.Y, self.Z, vmin=depth_extent[0], 
                                                       vmax=depth_extent[1],
                                                       cmap=cmap)
        else:
            plot = axes.imshow(self.Z, vmin=depth_extent[0], 
                                       vmax=depth_extent[1],
                                       extent=region_extent, 
                                       cmap=cmap,
                                       origin='lower')
        cbar = plt.colorbar(plot, ax=axes)
        cbar.set_label("Depth (m)")
        # levels = range(0,int(-numpy.min(Z)),500)

        # Plot coastlines
        if coastlines:
            axes.contour(self.X, self.Y, self.Z, levels=[0.0],colors='r')

        axes.set_xlim(region_extent[0:2])
        axes.set_ylim(region_extent[2:])

        return axes


    def project_unstructured(self, X_fill, Y_fill, Z_fill, extent=None,
                                   method='nearest', delta_limit=20.0, 
                                   no_data_value=-9999, buffer_length=100.0,
                                   proximity_radius=100.0, 
                                   resolution_limit=2000):
        r"""Project unstructured data on to regular grid.

        Function to project the unstructured data in the topo object onto a 
        structured grid.  Utilized a bounding box plus a buffer of size 
        *buffer_length* (meters) containing all data unless *extent* is not 
        None.  Then uses the fill data provided (*X_fill*, *Y_fill* and 
        *Z_fill*) to fill in the gaps in the unstructured data.  By default this
        is done by masking the fill data with the extents, the value 
        *no_data_value* and if *proximity_radius* (meters) is not 0, by a radius
        of *proximity_radius* from all grid points in the object.  Stores the 
        result in the *self.X*, *self.Y* and *self.Z* object attributes.  The
        resolution of the final grid is determined by calculating the minimum
        distance between all *x* and *y* data with a hard lower limit of 
        *delta_limit* (meters).

        :Input:
         - *extent* (tuple) - A tuple defining the rectangle of the sub-section.  
           Must be in the form (x lower,x upper,y lower, y upper).
         - *method* (string) - Method used for interpolation, valid methods are
           found in *scipy.interpolate.griddata*.  Default is *nearest*.
         - *delta_limit* (float) - Limit of finest horizontal resolution, 
           default is 20 meters.
         - *no_data_value* (float) - Value to use if no data was found to fill in a 
           missing value, ignored if `method = 'nearest'`. Default is `-9999`.
         - *buffer_length* (float) - Buffer around bounding box, only applicable
           when *extent* is None.  Default is *100.0* meters.
         - *proximity_radius* (float) - Radius every unstructured data point
           used to mask the fill data with.  Default is *100.0* meters.
         - *resolution_limit* (int) - Limit the number of grid points in a
           single dimension.  Raises a *ValueError* if the limit is violated.
           Default value is 

        """

        import scipy.interpolate as interpolate

        # Convert meter inputs to degrees
        mean_latitude = numpy.mean(self.y)
        buffer_degrees = dist_meters2latlong(buffer_length, 0.0, mean_latitude)[0]
        delta_degrees = dist_meters2latlong(delta_limit, 0.0, mean_latitude)[0]
        if proximity_radius > 0.0:
            proximity_radius_deg = dist_meters2latlong(proximity_radius, 0.0,
                                                        mean_latitude)[0]
            
        # Calculate new grid coordinates
        if extent is None:
            extent = [ numpy.min(self.x) - buffer_degrees, 
                       numpy.max(self.x) + buffer_degrees, 
                       numpy.min(self.y) - buffer_degrees, 
                       numpy.max(self.y) + buffer_degrees ]
        delta = max( min(numpy.min(numpy.abs(self.x[1:] - self.x[:-1])), 
                         numpy.min(numpy.abs(self.y[1:] - self.y[:-1])) ),
                    delta_degrees)
        N = ( numpy.ceil((extent[1] - extent[0]) / delta),
              numpy.ceil((extent[3] - extent[2]) / delta) )
        assert numpy.all(N[:] < numpy.ones((2)) * resolution_limit), \
               ValueError("Calculated resolution too high, N=%s!" % str(N))
        self._X, self._Y = numpy.meshgrid( 
                                     numpy.linspace(extent[0], extent[1], N[0]),
                                     numpy.linspace(extent[2], extent[3], N[1]))

        # Create extent mask
        extent_mask = extent[0] > X_fill
        extent_mask = numpy.logical_or(extent_mask,extent[1] < X_fill)
        extent_mask = numpy.logical_or(extent_mask,extent[2] > Y_fill)
        extent_mask = numpy.logical_or(extent_mask,extent[3] < Y_fill)
        
        # Create fill no-data value mask
        no_data_mask = numpy.logical_or(extent_mask, Z_fill == no_data_value)

        all_mask = numpy.logical_or(extent_mask, no_data_mask)

        # Create proximity mask
        if proximity_radius > 0.0:
        
            indices = (~all_mask).nonzero()
            for n in xrange(indices[0].shape[0]):
                i = indices[0][n]
                j = indices[1][n]
                all_mask[i,j] = numpy.any(numpy.sqrt((self.x - X_fill[i,j])**2 
                                                   + (self.y - Y_fill[i,j])**2)
                                             < proximity_radius_deg)

        X_fill_masked = numpy.ma.masked_where(all_mask, X_fill)
        Y_fill_masked = numpy.ma.masked_where(all_mask, Y_fill)
        Z_fill_masked = numpy.ma.masked_where(all_mask, Z_fill)    

        # Stick both the input data and fill data into arrays
        fill_points = numpy.column_stack((X_fill_masked.compressed(),
                                          Y_fill_masked.compressed()))
        points = numpy.concatenate((numpy.array([self.x, self.y]).transpose(), 
                                    fill_points))
        values = numpy.concatenate((self.z, Z_fill_masked.compressed()))

        # Use nearest-neighbor interpolation
        self._Z = interpolate.griddata(points, values, (self.X, self.Y), 
                                                                  method=method)

        self._extent = extent
        self._delta = delta
        self.unstructured = False


    def in_poly(self, polygon):
        r"""Mask points (x,y) that are not in the specified polygon.

        Uses simple ray casting algorithm for speed so beware of corner cases!

        Input
        -----
         - *polygon* (list) List of points that comprise the polygon.  Note that
           order of the points will effect if this works (positive versus negative
           winding order).  Points should be in counter-clockwise arrangement.

        Returns
        -------
         - *X_mask* (numpy.ma.MaskedArray) Masked array of X coordinates where those
           points outside of the polygon have been masked.
         - *Y* (numpy.ndarray) Coordinates in y direction in a meshgrid type of
           configuration.

        """
        raise NotImplemented("This function is not quite working yet, please "+\
                             "try again later")

        TOLERANCE = 1e-6

        # Flatten the input arrays to make this a bit easier
        x = self.X.flatten()
        y = self.Y.flatten()

        # Construct edges
        edges = []
        for edge in xrange(len(polygon) - 1):
            edges.append([polygon[edge], polygon[edge+1]])
        edges.append([polygon[-1], polygon[0]])

        # Check for intersections
        num_intersections = numpy.zeros(x.shape[0])

        for edge in edges:
            # Check for a vertical line
            if numpy.abs(edge[0][0] - edge[1][0]) < TOLERANCE:
                x_intersect = edge[0][0]        
            else:
                edge_slope = (edge[0][1] - edge[1][1]) / (edge[0][0] - edge[1][0])
                x_intersect = (y - edge[0][1]) / edge_slope + edge[0][0]

            num_intersections += (min(edge[0][1], edge[1][1]) <= y) * \
                                 (max(edge[0][1], edge[1][1]) >= y) * \
                                 (x_intersect <= x)
                                 

        # General intersection of two lines
        intersect = (numpy.mod(num_intersections, numpy.ones(x.shape) * 2) != 1)

        # Return masked arrays that are reshaped back to the input shapes
        return numpy.ma.masked_where(intersect, x, copy=False).reshape(X.shape), \
               numpy.ma.masked_where(intersect, y, copy=False).reshape(Y.shape)


    def replace_values(self, indices, value=numpy.nan, method='fill'):
        r"""Replace the values at *indices* by the specified method

        :Methods:
         - "fill"
         - "nearest"
        """

        # Average surrounding good data
        if method == 'fill':
            for index in indices:
                r = 0
                point_replaced = False
                while not point_replaced and r < max(self.Z.shape):
                    r = r + 1
                    i_range = range(max(0, index[0] - r), min(index[0] + r + 1, self.Z.shape[0]))
                    j_range = range(max(0, index[1] - r), min(index[1] + r + 1, self.Z.shape[1]))
                    num_points = 0
                    summation = 0.0
                    for i in i_range:
                        for j in j_range:
                            if (i,j) not in indices:
                                summation += self.Z[i,j]
                                num_points += 1
                    if num_points > 0:
                        Z[index[0], index[1]] = summation / num_points
                        point_replaced = True

        elif method == "nearest":
            pass

    def replace_no_data_values(self, method='fill'):
        r"""Replace *no_data_value* with other values as specified by *method*.

        self.no_data_value

        :Methods:
         - *fill* - Fill in all *no_data_value*s with *value*
         - *nearest* - Fill in *no_data_value*s with average of nearest
           neighbors.

        """
        raise NotImplemented("This functionality has not been added yet.")
        no_data_value_indices = (self.Z == self.no_data_value).nonzero()
        self.replace_values(no_data_value_indices, method=method)

        # nrows= shape(Z)[0]
        # ncols= shape(Z)[1]
        # npts = nrows*ncols

        # xi=X[0,:]
        # yi=Y[:,0]

        # X.np.reshape(npts)
        # Y.np.reshape(npts)
        # Z.np.reshape(npts)

        # ind=np.where(Z!=nodata_value)
        # X=X[ind]
        # Y=Y[ind]
        # Z=Z[ind]

        # ptsremove=npts-len(Z)
        # if ptsremove>0:
        #     print("Removing %s nodata_value points" % ptsremove)

        # Z = pylab.griddata(X,Y,Z,xi,yi)
        # (X,Y)=np.meshgrid(xi,yi)

        # griddata2topofile(X,Y,Z,outputfile,topotypeout,nodata_value,nodata_value)

    def smooth_data(self, indices, r=1):
        r"""Filter depth data at *indices* by averaging surrounding data.

        Surrounding data is considered within the ball of radius *r* in the 
        inf-norm.  Acts as a low-band pass filter and removes oscillatory data.

        :Input:
         - *indices* (list)
         - *r* (int) 

        :Output:
         None
        """

        index_range = [None, None]
        for index in indices:
            for n in xrange(2):
                index_range[n] = range(max(0, index[n] - r), 
                                       min(index[n] + r + 1, self.Z.shape[n]))
            num_points = 0
            summation = 0.0
            for i in index_range[0]:
                for j in index_range[1]:
                    summation += self.Z[i,j]
                    num_points += 1
            if num_points > 0:
                self.Z[index[0], index[1]] = summation / num_points

    def crop(self,filter_region):
        r"""Crop region to *filter_region*

        Create a new Topography object that is identical to this one but cropped
        to the region specified by filter_region

        :TODO:
         - Currently this does not work for unstructured data, could in principle
         - This could be a special case of in_poly although that routine could
           leave the resulting topography as unstructured effectively.
        """

        # Find indices of region
        region_index = [None, None, None, None]
        region_index[0] = (self.x >= filter_region[0]).nonzero()[0][0]
        region_index[1] = (self.x <= filter_region[1]).nonzero()[0][-1]
        region_index[2] = (self.y >= filter_region[2]).nonzero()[0][0]
        region_index[3] = (self.y <= filter_region[3]).nonzero()[0][-1]
        newtopo = Topography()

        newtopo._x = self._x[region_index[0]:region_index[1]]
        newtopo._y = self._y[region_index[2]:region_index[3]]
        if self._z is not None:
            newtopo._z = self._z[region_index[0]:region_index[1]]

        # Force regeneration of 2d coordinate arrays and extent if needed
        newtopo._X = None
        newtopo._Y = None
        newtopo._extent = None

        # Modify Z array as well
        newtopo._Z = self._Z[region_index[2]:region_index[3],
                          region_index[0]:region_index[1]]

        newtopo.unstructured = self.unstructured
        newtopo.topo_type = self.topo_type

        # print "Cropped to %s by %s array"  % (len(newtopo.x),len(newtopo.y))
        return newtopo
