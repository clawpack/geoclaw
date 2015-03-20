"""
Module for converting MOST data and output to GeoClaw format.

This is Work In Progress -- this module should be improved to use newer
topotools.

"""

import os, glob, re

def most2tt3(fname, fname2=None):
    """
    Converts MOST topo file to tt3 format.
    """
    f = open(fname).readlines()
    mn = f[0].split()
    ncols = int(mn[0])
    nrows = int(mn[1])

    xll = float(f[1]) 
    xup = float(f[ncols])
    #dx = float(f[2]) - xll
    dx = (xup - xll) / (ncols-1.)
    xll = xll - 360.

    yll = float(f[nrows+ncols])
    yup = float(f[ncols+1])
    #dy =  float(f[nrows+ncols-1]) - yll
    dy =  (yup - yll) / (nrows - 1.)
    if abs(dx-dy) > 1.e-6:
        print '*** WARNING: dx = ',dx,'  dy = ',dy
        yoffset = abs(yll + (nrows-1)*dx - yup) * 111e3
        print '*** This may cause offset in topography by  %s meters' % yoffset
        print '*** Consider using most2tt5'
    cellsize = dx

    # x = array([float(j) for j in f[1:ncols+1]])
    # y = array([float(j) for j in f[ncols+1:nrows+ncols+1]])

    if fname2 is None:
        fname2 = os.path.splitext(fname)[0] + '.tt3'

    f2 = open(fname2,'w')
    f2.write('%s ncols\n%s nrows\n%s xll\n%s yll\n%s cellsize\n99999 nodata_value\n' \
          % (ncols,nrows,xll,yll,cellsize))
    f2.writelines(f[nrows+ncols+1:])
    f2.close()
    print "Created ",fname2

def most2tt5(fname, fname2=None):
    """
    Converts MOST topo file to tt5 format.
    """
    from numpy import array, diff
    f = open(fname).readlines()
    mn = f[0].split()
    ncols = int(mn[0])
    nrows = int(mn[1])

    xll = float(f[1]) 
    xup = float(f[ncols])
    dx = (xup - xll) / (ncols-1.)
    xll = xll - 360.

    yll = float(f[nrows+ncols])
    yup = float(f[ncols+1])
    dy =  (yup - yll) / (nrows - 1.)
    print 'dx = ',dx,'  dy = ',dy

    x = array([float(j) for j in f[1:ncols+1]]) - 360.
    ddx = diff(x).max() - diff(x).min()
    if abs(ddx)/dx > 1e-3:
        print "*** Warning, dx not constant, ddx = ",ddx
        print "*** Data may not be equally spaced in x"
    y = array([float(j) for j in f[ncols+1:nrows+ncols+1]])
    ddy = diff(y).max() - diff(y).min()
    if abs(ddy)/dy > 1e-3:
        print "*** Warning, dy not constant, ddy = ",ddy
        print "*** Data may not be equally spaced in y"

    if fname2 is None:
        fname2 = os.path.splitext(fname)[0] + '.tt5'

    f2 = open(fname2,'w')
    f2.write('%s ncols\n%s nrows\n%s xll\n%s yll\n%s dx\n%s dy\n99999 nodata_value\n' \
          % (ncols,nrows,xll,yll,dx,dy))
    f2.writelines(f[nrows+ncols+1:])
    f2.close()
    print "Created ",fname2

def most2fortt(fnameprefix):
    """
    Converts MOST output files to fort.t files.
    """
    files = glob.glob(r'%s*' % fnameprefix)
    files.sort()
    s = r"%s(?P<hours>[0-9]*)h(?P<minutes>[0-9]*)m(?P<seconds>[0-9]*)s" \
          % fnameprefix
    regexp = re.compile(s)
    frameno = 1
    for fname in files:
        result = regexp.search(fname)
        try:
            hours = result.group("hours")
            minutes = result.group("minutes")
            seconds = result.group("seconds")
        except:
            print "*** Cannot parse fname: ",fname
            raise
        t = int(hours)*3600. + int(minutes)*60. + int(seconds)
        fortname = "fort.t" + str(frameno).zfill(4)
        f = open(fortname, 'w')
        f.write("%18.8e     time\n" % t)
        f.write("%5i                  meqn\n" % 1)
        f.write("%5i                  ngrids\n" % 1)
        f.write("%5i                  ndim\n" % 0)
        f.write("%5i                  maux\n" % 2)
        f.close()
        print "Created %s from %s at time t = %s" % (fortname, fname, t)
        frameno = frameno + 1


def most2fortq(fnameprefix):
    """
    Converts MOST output files to fort.q files.
    """
    files = glob.glob(r'%s*' % fnameprefix)
    files.sort()
    frameno = 1
    for fname in files:
        f = open(fname).readlines()
        mn = f[0].split()
        ncols = int(mn[0])
        nrows = int(mn[1])
        xll = float(f[1]) 
        dx = float(f[2]) - xll
        xll = xll - 360.
        yll = float(f[nrows+ncols])
        dy =  float(f[nrows+ncols-1]) - yll
        if abs(dx-dy) > 1.e-6:
            print '*** WARNING: dx = ',dx,'  dy = ',dy
        cellsize = dx
    
        fortname = 'fort.q' + str(frameno).zfill(4)
        f2 = open(fortname,'w')
        f2.write("%5i                  grid_number\n" % 1)
        f2.write("%5i                  AMR_level\n" % 1)
        f2.write("%5i                  mx\n" % ncols)
        f2.write("%5i                  my\n" % nrows)
        f2.write("%5i                  xlow\n" % xll)
        f2.write("%5i                  ylow\n" % yll)
        f2.write("%5i                  dx\n" % dx)
        f2.write("%5i                  dy\n" % dy)
        f2.write("\n")

        for k in range(len(f)-1, nrows+ncols, -1):
            for s in f[k].split():
                z = float(s)
                f2.write("%18.8e\n" % z)
        f2.close()
        print "Created %s from %s" % (fortname,fname)
        frameno += 1

def get_comMIT_topo(x1,x2,y1,y2, dx='1m', fname='most_topo.most'):
    r"""
    Get a topo file from the comMIT Bathymetry server, see
      http://sift.pmel.noaa.gov/ComMIT/bathymetry_help.html
    Inputs:
      - *x1, x2* (floats) extent of longitude (between 0 and 360)
      - *y1, y2* (floats) extent of latitude (between -90 and 90))
      - *dx* (float or str) cellsize in decimal degrees, 
             or string such as '1m', '30s' for arcminutes or arcseconds.
        Note that etopo1 bathymetry is used to construct this file and then
        linear interpolation is used if the resolution requested is finer
        than 1 minute.  But onshore topography comes from a different database
        and may be finer.
      - *fname* (str) File name to save as.

    The result is in MOST format.  
    You can convert to GeoClaw form using most2tt3, e.g.
        >>> most2geoclaw.most2tt3('most_topo.most')
    But note that this seems to create a file with topo_type == -3.

    Retrieving seems to only work for some values of x1,x2,y1,y2 and not others
    for reasons not yet understood!
      
    """
    import os, time
    import subprocess
    fname_most = os.path.splitext(fname)[0] + '.most'

    s = "wget -O " + fname_most + \
        " 'http://sift.pmel.noaa.gov/ComMIT/bathymetry/create.most" + \
        "?xmin=%s&xmax=%s&ymin=%s&ymax=%s&cellsize=%s'"  % (x1,x2,y1,y2,dx)

    wait_time = 3
    time.sleep(wait_time)
    subprocess.check_call(s,shell=True)
    lines = open(fname_most).readlines()
    attempts = 1
    while (len(lines) < 2) and (attempts < 5):
        print "Status: ",lines[0]
        print "Will retry in %s seconds..." % wait_time
        time.sleep(wait_time)
        subprocess.check_call(s,shell=True)
        lines = open(fname_most).readlines()
        attempts += 1
    if attempts<5:
        print "File downloaded with dimensions: ",lines[0]
    else:
        print "Tried %s times and file not ready... " % attempts
        print "     wait a bit and try same call again"



if __name__=='__main__':
    import sys
    most2tt3(sys.argv[1])
