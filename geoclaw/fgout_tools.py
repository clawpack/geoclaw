r"""
fgout_tools module: $CLAW/geoclaw/src/python/geoclaw/fgout_tools.py

Tools to specify and work with fgout grids, used to output GeoClaw solutions
on a fixed grid at a sequence of times, regardless of the AMR structure.

Includes:

- class FGoutFrame: used to hold a single frame of fgout output data
- class FGoutGrid: used to specify and store info about an fgout grid, with
                   methods to read and write info to fgout_grids.data
- function make_fgout_fcn_xy: Takes an FGoutFrame object and produces an
            interpolating function that can be evaluated for any (x,y).
- function make_fgout_fcn_xyt: Takes 2 FGoutFrame objects and produces an
            interpolating function that can be evaluated for any (x,y,t)
            at intermediate times.
- function write_netcdf: Write a specified set of qoi's from a list of
            fgout frames, as a single netCDF file
- function read_netcdf: Read a netCDF file and return a list of fgout frames,
            assuming the file contains all the qoi's needed to reconstruct q.
- function read_netcdf_arrays: Read a netCDF file and extract the 
            requested quantities of interest as numpy arrays.
- print_netcdf_info: Print info about the contents of a netCDF file containing
            some fgout frames.
"""

import os
from numpy import sqrt, ma, mod
import numpy
from six.moves import range


class FGoutFrame(object):

    """
    Class to hold a single frame of fgout data at one output time.
    Several attributes are defined as properties that can be evaluated
    and stored only when needed by the user.
    """

    def __init__(self, fgout_grid, frameno=None):
        self.fgout_grid = fgout_grid
        self.frameno = frameno
        self.t = None

        # private attributes for those that are only created if
        # needed by the user:
        self._h = None
        self._hu = None
        self._hv = None
        self._eta = None
        self._B = None
        self._u = None
        self._v = None
        self._s = None
        self._hss = None

    # Define shortcuts to attributes of self.fgout_grid that are the same
    # for all frames (e.g. X,Y) to avoid storing grid for every frame.
    
    @property
    def x(self):
        return self.fgout_grid.x

    @property
    def y(self):
        return self.fgout_grid.y
        
    @property
    def X(self):
        return self.fgout_grid.X
        
    @property
    def Y(self):
        return self.fgout_grid.Y
        
    @property
    def delta(self):
        return self.fgout_grid.delta

    @property
    def extent_centers(self):
        return self.fgout_grid.extent_centers
        
    @property
    def extent_edges(self):
        return self.fgout_grid.extent_edges
        
    @property
    def drytol(self):
        return self.fgout_grid.drytol
            
    # Define attributes such as h as @properties with lazy evaluation:
    # the corresponding array is created and stored only when first
    # accessed by the user.  Those not needed are not created.

    @property
    def h(self):
        """depth"""
        if self._h is None:
            #print('+++ setting _h...')
            self._h = self.q[0,:,:]
        #print('+++ getting _h...')
        return self._h

    @property
    def hu(self):
        """momentum h*u"""
        if self._hu is None:
            self._hu = self.q[1,:,:]
        return self._hu

    @property
    def u(self):
        """speed u, computed as hu/h or set to 0 if h<self.drytol"""
        if self._u is None:
            self._u = numpy.divide(self.hu, self.h,\
                              out=numpy.zeros(self.h.shape, dtype=float), \
                              where=(self.h>self.drytol))
        return self._u

    @property
    def hv(self):
        """momentum h*v"""
        if self._hv is None:
            self._hv = self.q[2,:,:]
        return self._hv

    @property
    def v(self):
        """speed v, computed as hv/h or set to 0 if h<self.drytol"""
        if self._v is None:
            self._v = numpy.divide(self.hv, self.h,\
                              out=numpy.zeros(self.h.shape, dtype=float), \
                              where=(self.h>self.drytol))
        return self._v

    @property
    def eta(self):
        """surface eta = h+B"""
        if self._eta is None:
            self._eta = self.q[-1,:,:]
        return self._eta

    @property
    def B(self):
        """topography"""
        if self._B is None:
            self._B = self.q[-1,:,:] - self.q[0,:,:]
        return self._B

    @property
    def s(self):
        """speed s = sqrt(u**2 + v**2)"""
        if self._s is None:
            self._s = numpy.sqrt(self.u**2 + self.v**2)
        return self._s

    @property
    def hss(self):
        """momentum flux h*s**2"""
        if self._hss is None:
            self._hss = self.h * self.s**2
        return self._hss


class FGoutGrid(object):

    """
    New class introduced in 5.9.0 to keep store information both about the
    fgout input data and the output generated by a GeoClaw run.
    """

    def __init__(self,fgno=None,outdir=None,output_format=None):


        # GeoClaw input values:
        self.id = ''  # identifier, optional
        self.point_style = 2  # only option currently supported
        self.npts = None
        self.nx = None
        self.ny = None
        self.tstart =  None
        self.tend = None
        self.nout = None
        self.fgno = fgno
        self.outdir = outdir
        self.output_format = output_format
        
        self.drytol = 1e-3          # used for computing u,v from hu,hv

        # private attributes for those that are only created if
        # needed by the user:
        
        self._X = None
        self._Y = None
        self._x = None
        self._y = None
        self._delta = None           # (dx,dy)
        self._extent_centers = None  # defined by fgout points
        self._extent_edges = None    # extended so points are cell centers

        self._plotdata = None

    @property
    def x(self):
        """1D array x of longitudes (cell centers)"""
        if self._x is None:
            dx = (self.x2 - self.x1)/self.nx
            self._x = numpy.linspace(self.x1+dx/2, self.x2-dx/2, self.nx)
        return self._x

    @property
    def y(self):
        """1D array y of latitudes (cell centers)"""
        if self._y is None:
            dy = (self.y2 - self.y1)/self.ny
            self._y = numpy.linspace(self.y1+dy/2, self.y2-dy/2, self.ny)
        return self._y
        
    @property
    def X(self):
        """2D array X of longitudes (cell centers)"""
        if self._X is None:
            self._X = numpy.meshgrid(self.x, self.y)[0].T
        return self._X

    @property
    def Y(self):
        """2D array Y of latitudes (cell centers)"""
        if self._Y is None:
            self._Y = numpy.meshgrid(self.x, self.y)[1].T
        return self._Y

    @property
    def delta(self):
        if self._delta is None:
            dx = (self.x2 - self.x1)/self.nx
            dy = (self.y2 - self.y1)/self.ny
            self._delta = (dx,dy)
        return self._delta
        
    @property
    def extent_centers(self):
        """Extent of cell centers [xmin,xmax,ymin,ymax]"""
        if self._extent_centers is None:
            self._extent_centers = [self.x.min(), self.x.max(),\
                                    self.y.min(), self.y.max()]
        return self._extent_centers
        
    @property
    def extent_edges(self):
        """Extent of cell edges [xmin,xmax,ymin,ymax]"""
        if self._extent_edges is None:
            dx,dy = self.delta
            self._extent_edges = [self.x.min()-dx/2, self.x.max()+dx/2,\
                                    self.y.min()-dy/2, self.y.max()+dy/2]
        return self._extent_edges
        

    # Create plotdata of class clawpack.visclaw.ClawPlotData
    # only when needed for reading GeoClaw output, 
    # since this must be done after fgno, outdir, output_format
    # have been specified:

    @property
    def plotdata(self):
        if self._plotdata is None:
            self._plotdata = self.set_plotdata()
        return self._plotdata

    def set_plotdata(self):
        """
        Create a plotdata, assuming attributes fgno, outdir, output_format
        have all been set.
        """
        from clawpack.visclaw.data import ClawPlotData

        assert self.fgno is not None, '*** fgno must be set'
        assert self.outdir is not None, '*** outdir must be set'
        assert self.output_format is not None, '*** output_format must be set'
        if 0:
            print('+++ creating plotdata for fgno=%i, outdir=%s, format=%s' \
                  % (self.fgno,self.outdir,self.output_format))

        plotdata = ClawPlotData()
        plotdata.outdir = self.outdir
        plotdata.format = self.output_format
        plotdata.file_prefix = 'fgout%s' % str(self.fgno).zfill(4)
        if self.output_format[:6]=='binary':
            # could be 'binary', 'binary32' or 'binary64'
            file_prefix_str = plotdata.file_prefix + '.b'
        else:
            file_prefix_str = plotdata.file_prefix + '.q'
        
        # Contrary to usual default, set save_frames to False so that all
        # the fgout frames read in are not saved in plotdata.framesoln_dict.
        # Otherwise this might be a memory hog when making an animation with
        # many frames on a fine fgout grid.
        plotdata.save_frames = False

        return plotdata


    def read_fgout_grids_data(self, fgno=None, data_file='fgout_grids.data'):
        """
        Read input info for fgout grid number fgno from the data file
        fgout_grids.data, which should have been created by setrun.py.
        This file now contains info about all fgout grids.
        """
        if fgno is not None:
            self.fgno = fgno
        assert self.fgno is not None, '*** fgno must be set'

        with open(data_file) as filep:
            lines = filep.readlines()
        fgout_input = None
        for lineno,line in enumerate(lines):
            if 'fgno' in line:
                if int(line.split()[0]) == self.fgno:
                    fgout_input = lines[lineno+1:]
                    #print('Found line %i: %s' % (lineno,line))
                    break

        if fgout_input is None:
            raise ValueError('fgout grid fgno = %i not found in %s' \
                             % (fgno, data_file))

        self.tstart = float(fgout_input[0].split()[0])
        self.tend = float(fgout_input[1].split()[0])
        self.nout = int(fgout_input[2].split()[0])
        self.point_style = point_style = int(fgout_input[3].split()[0])
        output_format = int(fgout_input[4].split()[0])
        if output_format == 1:
            self.output_format = 'ascii'
        elif output_format == 3:
            self.output_format = 'binary'
        print('Reading input for fgno=%i, point_style = %i ' \
                % (self.fgno, self.point_style))
        if point_style == 2:
            self.nx = nx = int(fgout_input[5].split()[0])
            self.ny = ny = int(fgout_input[5].split()[1])
            self.x1 = float(fgout_input[6].split()[0])
            self.y1 = float(fgout_input[6].split()[1])
            self.x2 = float(fgout_input[7].split()[0])
            self.y2 = float(fgout_input[7].split()[1])
        else:
            raise NotImplementedError("fgout not implemented for point_style %i" \
                % point_style)


    def write_to_fgout_data(self, fid):
        """
        Convert fgout data specified in setrun.py to file `fgout_grids.data`
        read in by GeoClaw fortran code.
        """

        print("\n---------------------------------------------- ")
        assert self.point_style is not None, 'Need to set point_style'

        point_style = self.point_style
        if point_style not in [2]:
            errmsg = 'point_style %s is not implement, only point_style==2' \
                                    % point_style
            raise NotImplementedError(errmsg)
                
        if self.output_format == 'ascii':
            output_format = 1
        elif self.output_format == 'binary32':
            output_format = 2
        elif self.output_format in ['binary','binary64']:
            output_format = 3
        else:
            errmsg = "fgout output_format must be ascii, binary32, or binary64"
            raise NotImplementedError(errmsg)

        assert self.tstart is not None, 'Need to set tstart'
        assert self.tend is not None, 'Need to set tend'
        assert self.nout is not None, 'Need to set nout'
        assert self.point_style is not None, 'Need to set point_style'
        
        # write header, independent of point_style:
        #fid = open(self.input_file_name,'w')
        fid.write("\n")
        fid.write("%i                           # fgno\n" % self.fgno)
        fid.write("%16.10e            # tstart\n"  % self.tstart)
        fid.write("%16.10e            # tend\n"  % self.tend)
        fid.write("%i %s           # nout\n" % (self.nout, 11*" "))
        fid.write("%i %s              # point_style\n" \
                            % (self.point_style,12*" "))
        fid.write("%i %s              # output_format\n" \
                            % (output_format,12*" "))

        print('fgout grid %i has point_style = %i' % (self.fgno, point_style))


        if point_style == 2:
            # 2d grid of points
            x1,x2 = self.x1, self.x2
            y1,y2 = self.y1, self.y2
            nx,ny = self.nx, self.ny

            dx = (x2-x1)/nx   # x1,x2 are cell edges
            dy = (y2-y1)/ny   # y1,y2 are cell edges

            npts = nx*ny

            fid.write("%i  %i %s          # nx,ny\n" \
                                % (nx,ny,10*" "))
            fid.write("%16.10e   %20.10e            # x1, y1\n" % (x1,y1))
            fid.write("%16.10e   %20.10e            # x2, y2\n" % (x2,y2))

            print("   specifying fgout grid with shape %i by %i, with  %i points" \
                    % (nx,ny,npts))
            print("   lower left  = (%15.10f,%15.10f)" % (x1,y1))
            print("   upper right = (%15.10f,%15.10f)" % (x2,y2))
            print("   dx = %15.10e,  dy = %15.10e" % (dx,dy))
            


    def read_frame(self, frameno):
        """
        Read a single frame of fgout data.
        """

        from datetime import timedelta

        try:
            fr = self.plotdata.getframe(frameno)
        except:
            print('*** Could not read fgout grid %i frame %i from %s' \
                 % (self.fgno,frameno,self.plotdata.outdir))
            raise
        state = fr.states[0]  # only 1 AMR grid
        patch = state.patch

        fgout_frame = FGoutFrame(self, frameno)
        fgout_frame.fgout_grid = self

        fgout_frame.q = state.q

        fgout_frame.t = state.t
    
        fgout_frame.frameno = frameno
        
        X,Y = patch.grid.p_centers[:2]
        
        try:
            fgoutX = self.X
            fgoutY = self.Y
        except:
            # presumably x,y, etc. were not set for this fgout_grid
            # (e.g. by read_fgout_grids_data) so infer from this frame and set
            
            # reset most attributes including private _x etc to None:
            self.__init__(fgno=self.fgno,outdir=self.outdir,
                          output_format=self.output_format)
                          
            print('    Setting grid attributes based on Frame %i:' % frameno)
            x = X[:,0]
            y = Y[0,:]
            dx = x[1] - x[0]
            dy = y[1] - y[0]
            self.x1 = x[0] - dx/2
            self.x2 = x[-1] + dx/2
            self.y1 = y[0] - dy/2
            self.y2 = y[-1] + dy/2
            self.nx = len(x)
            self.ny = len(y)
            fgoutX = self.X  # will be computed from info above
            fgoutY = self.Y  # will be computed from info above
            print('      Grid edges: ',self.extent_edges)
            print('      with %i cells in x and %i cells in y' \
                    % (self.nx,self.ny))
            
            
        if not numpy.allclose(X, fgoutX):
            errmsg = '*** X read from output does not match fgout_grid.X'
            raise ValueError(errmsg)

        if not numpy.allclose(Y, fgoutY):
            errmsg = '*** Y read from output does not match fgout_grid.Y'
            raise ValueError(errmsg)
                
        return fgout_frame


# ========================
# Functions for interpolating from fgout grid to arbitrary points,
# useful for example if using velocity field to model particle/debris motion

def make_fgout_fcn_xy(fgout, qoi, method='nearest',
                       bounds_error=False, fill_value=numpy.nan):
    """
    Create a function that can be called at (x,y) and return the qoi
    interpolated in space from the fgout array.
    
    qoi should be a string (e.g. 'h', 'u' or 'v') corresponding to 
    an attribute of fgout.
    
    The function returned takes arguments x,y that can be floats or 
    (equal length) 1D arrays of values that lie within the spatial
    extent of fgout. 
    
    bounds_error and fill_value determine the behavior if (x,y) is not in 
    the bounds of the data, as in scipy.interpolate.RegularGridInterpolator.
    """
    
    from scipy.interpolate import RegularGridInterpolator
    
    try:
        q = getattr(fgout,qoi)
    except:
        print('*** fgout missing attribute qoi = %s?' % qoi)
        
    err_msg = '*** q must have same shape as fgout.X\n' \
            + 'fgout.X.shape = %s,   q.shape = %s' % (fgout.X.shape,q.shape)
    assert fgout.X.shape == q.shape, err_msg
    
    x1 = fgout.X[:,0]
    y1 = fgout.Y[0,:]
    fgout_fcn1 = RegularGridInterpolator((x1,y1), q, method=method,
                bounds_error=bounds_error, fill_value=fill_value)

    def fgout_fcn(x,y):
        """
        Function that can be evaluated at single point or arrays (x,y).
        """
        from numpy import array, vstack
        xa = array(x)
        ya = array(y)
        xyout = vstack((xa,ya)).T
        qout = fgout_fcn1(xyout)
        if len(qout) == 1:
            qout = qout[0]  # return scalar
        return qout

    return fgout_fcn


def make_fgout_fcn_xyt(fgout1, fgout2, qoi, method_xy='nearest',
                       method_t='linear', bounds_error=False,
                       fill_value=numpy.nan):
    """
    Create a function that can be called at (x,y,t) and return the qoi
    interpolated in space and time between the two frames fgout1 and fgout2.
    
    qoi should be a string (e.g. 'h', 'u' or 'v') corresponding to 
    an attribute of fgout.
    
    method_xy is the method used in creating the spatial interpolator,
    and is passed to make_fgout_fcn_xy.
    
    method_t is the method used for interpolation in time, currently only
    'linear' is supported, which linearly interpolates.
    
    bounds_error and fill_value determine the behavior if (x,y,t) is not in 
    the bounds of the data.
    
    The function returned takes arguments x,y (floats or equal-length 1D arrays)
    of values that lie within the spatial extent of fgout1, fgout2
    (which are assumed to cover the same uniform grid at different times)
    and t should be a float that lies between fgout1.t and fgout2.t. 
    """
    
    assert numpy.allclose(fgout1.X, fgout2.X), \
                            '*** fgout1 and fgout2 must have same X'
    assert numpy.allclose(fgout1.Y, fgout2.Y), \
                            '*** fgout1 and fgout2 must have same Y'
    
    t1 = fgout1.t
    t2 = fgout2.t
    #assert t1 < t2, '*** expected fgout1.t < fgout2.t'
    
    fgout1_fcn_xy = make_fgout_fcn_xy(fgout1, qoi, method=method_xy,
                       bounds_error=bounds_error, fill_value=fill_value)
    fgout2_fcn_xy = make_fgout_fcn_xy(fgout2, qoi, method=method_xy,
                       bounds_error=bounds_error, fill_value=fill_value)
                
    def fgout_fcn(x,y,t):
        """
        Function that can be evaluated at single point or arrays (x,y)
        at a single time t.
        """
        from numpy import array, ones
        xa = array(x)
        ya = array(y)
        tol = 1e-6  # to make sure it works ok when called with t=t1 or t=t2
        if t1-tol <= t <= t2+tol:
            alpha = (t-t1)/(t2-t1)
        elif bounds_error:
            errmsg = '*** argument t=%g should be between t1=%g and t2=%g' \
                     % (t,t1,t2)
            raise ValueError(errmsg)
        else:
            qout = fill_value * ones(xa.shape)
            return qout

        qout1 = fgout1_fcn_xy(x,y)
        qout2 = fgout2_fcn_xy(x,y)

        if method_t == 'linear':
            if t1 <= t <= t2:
                alpha = (t-t1)/(t2-t1)
                qout = (1-alpha)*qout1 + alpha*qout2
        else:
            raise NotImplementedError('method_t = %s not supported' % method_t)
            
        return qout
        
    return fgout_fcn

# ===============================
# Functions for writing a set of fgout frames as a netCDF file, and
# reading such a file:    
    
def write_netcdf(fgout_frames, fname_nc='fgout_frames.nc',
                 qois = ['h','hu','hv','eta'], datatype='f4',
                 include_B0=False, include_Bfinal=False,
                 description='', verbose=True):
    """
    Write a list of fgout frames (at different times on the same rectangular
    grid) to a single netCDF file, with some metadata and the topography,
    if desired.
    
    fgout_frames should be a list of FGoutFrame objects, all of the same size
    and at increasing times.
    
    fname_nc is the name of the file to write.
    
    qois is a list of strings, the quantities of interest to include in file.
        This could include any of: 

            'h', 'eta', 'hu', 'hv', 'u', 'v', 's', 'hss', 'B'.

        All other quantities can be computed from h, hu, hv, eta,
        the original fgout variables from GeoClaw, but for some applications
        you might only want to save 'h' and 's', for example.
    
    datatype should be 'f4' [default] or 'f8', specifying bytes per qoi value.
        'f8' has full precision of the original data, but the file will be
        twice as large and may not be needed for downstream applications.
        
    Note that the topography B = eta - h, so it is not necessary to store all
    three of these.  Also, B is often the same for all frames, so rather than
    storing B at each frame as a qoi, two other options are also provided
    (and then storing eta or h for all frames allows calculating the other):
    
    include_Bfinal: If True, include the topography B array from the final frame
    as the Bfinal array.
    
    include_B0: If True, include the topography B array from the first frame
    as the B0 array.  This is only useful if, e.g., the first frame is initial 
    topography before co-seismic deformation, and at later times the topography
    is always equal to Bfinal. 
    
    `description` is a string that will be added as metadata.
    A metadata field `history` will also be added, which includes the 
    time the file was created and the path to the directory where it was made. 
    """    
    
    import netCDF4
    from datetime import datetime, timedelta
    from cftime import num2date, date2num
    import time
    timestr = time.ctime(time.time())  # current time for metadata
    
    fg_times = numpy.array([fg.t for fg in fgout_frames])
    
    if verbose:
        print('Creating %s with fgout frames at times: ' % fname_nc)
        print(fg_times)
    
    fg0 = fgout_frames[0]
    x = fg0.x
    y = fg0.y
    
    xs = numpy.array([fg.x for fg in fgout_frames])
    ys = numpy.array([fg.y for fg in fgout_frames])
    # assert same for all times
    
    units = {'h':'meters', 'eta':'meters', 'hu':'m^2/s', 'hv':'m^2/s',
             'u':'m/s', 'v':'m/s', 's':'m/s', 'hss':'m^3/s^2', 'B':'meters'}

    with netCDF4.Dataset(fname_nc, 'w') as rootgrp:

        rootgrp.description = description
        rootgrp.history = "Created " + timestr
        rootgrp.history += " in %s;  " % os.getcwd()

        lon = rootgrp.createDimension('lon', len(x))
        longitudes = rootgrp.createVariable('lon','f8',('lon',))
        longitudes[:] = x
        longitudes.units = 'degrees_east'

        lat = rootgrp.createDimension('lat', len(y))
        latitudes = rootgrp.createVariable('lat','f8',('lat',))
        latitudes[:] = y
        latitudes.units = 'degrees_north'
        
        time = rootgrp.createDimension('time', len(fg_times))
        times = rootgrp.createVariable('time','f8',('time',))
        times[:] = fg_times
        times.units = 'seconds'
        
        if 0:
            # Could make times be datetimes relative to some event time, e.g.:
            times.units = 'seconds since 1700-01-26 21:00:00.0'
            times.calendar = 'gregorian'
            dates = [datetime(1700,1,26,21) + timedelta(seconds=ss) \
                        for ss in fg_times]
            times[:] = date2num(dates,units=times.units,calendar=times.calendar)
        
        if include_B0:
            B0 = rootgrp.createVariable('B0',datatype,('lon','lat',))
            B0[:,:] = fg0.B
            B0.units = 'meters'

        if include_Bfinal:
            fg_final = fgout_frames[-1]
            Bfinal = rootgrp.createVariable('Bfinal',datatype,('lon','lat',))
            Bfinal[:,:] = fg_final.B
            Bfinal.units = 'meters'
        
        for qoi in qois:
            qoi_frames = [getattr(fgout,qoi) for fgout in fgout_frames]
            qoi_var = rootgrp.createVariable(qoi,datatype,('time','lon','lat',))
            qoi_var[:,:,:] = qoi_frames
            qoi_var.units = units[qoi]

def get_as_array(var, rootgrp, verbose=True):
    """
    Utility function to retrieve variable from netCDF file and convert to 
    numpy array.
    """
    a = rootgrp.variables.get(var, None)
    if a is not None:
        if verbose: print('    Loaded %s with shape %s' % (var,repr(a.shape)))
        return numpy.array(a)
    else:
        if verbose: print('    Did not find %s' % var)
        return None    


def read_netcdf_arrays(fname_nc, qois, verbose=True):
    """
    Read a netCDF file and extract the quantities of interest denoted by
    strings in the list qois, which can include:

        'h', 'eta', 'hu', 'hv', 'u', 'v', 's', 'hss', 'B'.

    qois can also include the time-independent 'B0' and/or 'Bfinal'
    
    Returns 

        x, y, t, qoi_arrays

    where x,y define the longitude, latitudes, t is the times of the frames,
    and qoi_arrays is a dictionary indexed by the strings from qois.
    
    Each dict element is an array with shape (len(t), len(x), len(y))
    for time-dependent qoi's, or (len(x), len(y)) for B0 or Bfinal,
    or None if that qoi was not found in the netCDF file.
    """
     
    import netCDF4

    with netCDF4.Dataset(fname_nc, 'r') as rootgrp:
        if verbose:
            print('Reading data to fgout frames from nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('History:  ', rootgrp.history)

        x = get_as_array('lon', rootgrp, verbose)
        y = get_as_array('lat', rootgrp, verbose)
        t = get_as_array('time', rootgrp, verbose)
        
        vars = list(rootgrp.variables)
        
        qoi_arrays = {}
        for qoi in qois:
            qoi_array = get_as_array(qoi, rootgrp, verbose)
            qoi_arrays[qoi] = qoi_array
            
    return x,y,t,qoi_arrays
            
            

def read_netcdf(fname_nc, fgout_grid=None, verbose=True):
    """
    Read a netCDF file and return a list of FGoutFrame instances.
    This will only be possible if the netCDF file contains at least
    the qoi's 'h','hu','hv','eta' required to reconstruct the q array
    as output by GeoClaw.
    """
    
    import netCDF4

    with netCDF4.Dataset(fname_nc, 'r') as rootgrp:
        if 1:
            print('Reading data to fgout frames from nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('History:  ', rootgrp.history)

        x = get_as_array('lon', rootgrp, verbose)
        y = get_as_array('lat', rootgrp, verbose)
        t = get_as_array('time', rootgrp, verbose)

        vars = list(rootgrp.variables)
        

        for qoi in ['h','hu','hv','eta']:
            errmsg = '*** Cannot reconstruct fgout frame without %s' % qoi
            assert qoi in vars, errmsg
        h = get_as_array('h', rootgrp, verbose)
        hu = get_as_array('hu', rootgrp, verbose)
        hv = get_as_array('hv', rootgrp, verbose)
        eta = get_as_array('eta', rootgrp, verbose)
        
        if (x is None) or (y is None):
            print('*** Could not create grid')
        else:
            X,Y = numpy.meshgrid(x,y)
        
        fgout_frames = []
                
        for k in range(eta.shape[0]):
            fgout = FGoutFrame(fgout_grid=fgout_grid, frameno=k)
            fgout.x = x
            fgout.y = y
            fgout.t = t[k]
            fgout.q = numpy.empty((4,eta.shape[1],eta.shape[2]))
            fgout.q[0,:,:] = h[k,:,:]
            fgout.q[1,:,:] = hu[k,:,:]
            fgout.q[2,:,:] = hv[k,:,:]
            fgout.q[3,:,:] = eta[k,:,:]
            fgout.X = X
            fgout.Y = Y
            fgout_frames.append(fgout)
            
        print('Created fgout_frames as list of length %i' % len(fgout_frames))
        
    return fgout_frames

def print_netcdf_info(fname_nc):
    """
    Print out info about the contents of a netCDF file contining fgout frames,
    written using write_netcdf.
    """
    import netCDF4

    with netCDF4.Dataset(fname_nc, 'r') as rootgrp:
        x = get_as_array('lon', rootgrp, verbose=False)
        y = get_as_array('lat', rootgrp, verbose=False)
        t = get_as_array('time', rootgrp, verbose=False)

        vars = list(rootgrp.variables)
        
        print('===================================================')
        print('netCDF file %s contains:' % fname_nc)
        print('description: \n', rootgrp.description)
        print('history: \n', rootgrp.history)
        print('%i longitudes from %.6f to %.6f' % (len(x),x[0],x[-1]))
        print('%i latitudes from %.6f to %.6f' % (len(y),y[0],y[-1]))
        print('%i times from %.3f to %.3f' % (len(t),t[0],t[-1]))
        print('variables: ',vars)
        print('===================================================')
    
        
