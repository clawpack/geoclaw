r"""
kmltools module: $CLAW/geoclaw/src/python/geoclaw/kmltools.py

Tools to make kml files to overlay on Google Earth.
Note that color is in KML format, BGR with 2 hex digits for each, e.g.

  FF0000 is blue, 00FF00 is green,  0000FF is red, 00FF00 is yellow.

Actually it's an 8 hex digit number, where the first two digits are
transparency, but in this module these default to 'FF' (but you can specify
the full 8 digits if you want it transparent).

:Functions:
 - deg2dms - convert decimal degrees to (degrees, minutes, seconds)
 - regions2kml - create a kml outline for each regions specified in setrun
 - box2kml - create a kml outline from a rectangular box
 - quad2kml - create a kml outline for an arbitrary quadrilateral
 - poly2kml - create a kml outline for an arbitrary polygon
 - line2kml - create a kml line connecting 2 points
 - gauges2kml - create a kml marker for each gauge specified in setrun
 - topo2kml - create a kml outline for each topo grid specified in setrun
 - dtopo2kml - create a kml outline for each dtopo grid specified in setrun
 - fgmax2kml - create a kml outline for each fgmax grid specified in setrun
 - fgout2kml - create a kml outline for each fgout grid specified in setrun
 - make_input_data_kmls - make kml files for many things specified in setrun
 - pcolorcells_for_kml - version of pcolormesh with appropriate dpi and size
 - png2kml - create kml file wrapping a png figure to be viewed on GE
 - kml_build_colorbar - create a colorbar to display on GE
 - topo2kmz - create kmz file showing onshore and offshore topography
 - transect2kml - create kml file showing a set of points on a transect
 - dtopo_contours2kmz - create a kmz file containing contour plots of dtopo(s)
 - kml_header - used internally
 - kml_footer - used internally
 - kml_region - used internally
 - kml_gauge - used internally
 - kml_png - used internally
 - kml_cb - used internally
"""

try:
    from importlib import reload
except:
    pass # assume python2, already has reload

import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
import os

def f2s(x, num_digits=6):
    r"""
    Convert float to string in fixed point notation with at most
    *num_digits* digits of precision and trailing zeros removed,
    for printing nicely in kml description boxes.
    """
    format = '%' + '.%sf' % num_digits
    s = (format % x).rstrip('0')
    return s

def deg2dms(dy):
    r"""
    Convert decimal degrees to tuple (degrees, minutes, seconds)
    """

    from numpy import floor
    dy_deg = floor(dy)
    dy_min = floor((dy-dy_deg)*60.)
    dy_sec = (dy-dy_deg-dy_min/60.)*3600.
    return dy_deg,dy_min,dy_sec


def regions2kml(rundata=None,fname='regions.kml',verbose=True,combined=True):

    """
    Create a KML box for each AMR region specified for a GeoClaw run.

    :Inputs:

      - *rundata* - an object of class *ClawRunData* or None

        If *rundata==None*, try to create based on executing function *setrun*
        from the `setrun.py` file in the current directory.

      - *fname* (str) - resulting kml file.

      - *verbose* (bool) - If *True*, print out info about each region found

      - *combined* (bool) - If *True*, combine into single kml file with
        name given by *fname*.  This is the default.
        If False, *fname* is ignored and individual files are created for
        each region with names are Domain.kml, Region00.kml, etc.
        These will show up separately in GoogleEarth so they can be turned
        on or off individually.

    First create a box for the entire domain (in red) and then a box
    for each region (in white).

    :Example:

        >>> from clawpack.geoclaw import kmltools
        >>> kmltools.regions2kml()

    is equivalent to:

        >>> from clawpack.geoclaw import kmltools
        >>> from setrun import setrun
        >>> rundata = setrun()
        >>> kmltools.regions2kml(rundata)

    By default this creates a file named *regions.kml* that can be opened in
    Google Earth.

    """

    from numpy import cos,pi,floor

    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    fname_combined = 'FlagRegions.kml'

    clawdata = rundata.clawdata
    x1,y1 = clawdata.lower[0:]
    x2,y2 = clawdata.upper[0:]
    description = "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s\n" % (f2s(y1),f2s(y2))

    mx,my = clawdata.num_cells[0:]
    dx = (x2-x1)/float(mx)
    dx_meters = dx*111e3*cos(pi*0.5*(y1+y2)/180.)
    dy = (y2-y1)/float(my)
    dy_meters = dy*111e3
    if verbose:
        print("Domain:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))
    dx_deg,dx_min,dx_sec = deg2dms(dx)
    dy_deg,dy_min,dy_sec = deg2dms(dy)
    #print "Level 1 resolution:  dx = %g deg, %g min, %g sec = %g meters" \
    #   % (dx_deg,dx_min,dx_sec,dx_meters)
    levtext = "Level 1 resolution:  dy = %g deg, %g min, %g sec = %g meters\n" \
        % (dy_deg,dy_min,dy_sec,dy_meters)
    if verbose:
        print(levtext)
    description = description + levtext

    amr_levels_max = rundata.amrdata.amr_levels_max
    refinement_ratios_y = rundata.amrdata.refinement_ratios_y
    num_ref_ratios = len(refinement_ratios_y)
    if amr_levels_max > num_ref_ratios+1:
        raise IOError("*** Too few refinement ratios specified for " \
            + "amr_levels_max = %i" % amr_levels_max)
    dy_levels = (num_ref_ratios+1) * [dy]
    for k,r in enumerate(refinement_ratios_y):
        level = k+2
        dy = dy_levels[k] / r
        dy_levels[k+1] = dy
        dy_meters = dy*111e3
        dy_deg,dy_min,dy_sec = deg2dms(dy)
        levtext = "Level %s resolution:  dy = %g deg, %g min, %g sec = %g meters  (refined by %i)\n" \
                % (level,dy_deg,dy_min,dy_sec,dy_meters,r)
        if verbose:
            print(levtext)
        description = description + levtext

    if verbose:
        print("Allowing maximum of %i levels" % amr_levels_max)

    elev = 0.
    if not combined:
        fname = 'Domain.kml'
    else:
        fname = fname_combined

    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = 'Computational Domain'
    mapping['desc'] = description
    mapping['color'] = "0000FF"  # red
    mapping['width'] = 2

    region_text = kml_region(mapping)
    kml_text = kml_text + region_text

    if not combined:
        kml_text = kml_text + kml_footer()
        kml_file = open(fname,'w')
        kml_file.write(kml_text)
        kml_file.close()
        if verbose:
            print("Created ",fname)



    regions = rundata.regiondata.regions
    if len(regions)==0 and verbose:
        print("No regions found in setrun.py")


    for rnum,region in enumerate(regions):
        if not combined:
            fname = 'Region_%s.kml' % str(rnum).zfill(2)
            kml_text = kml_header(fname)

        minlevel,maxlevel = region[0:2]
        t1,t2 = region[2:4]
        x1,x2,y1,y2 = region[4:]

        if verbose:
            print("Region %i: %10.6f  %10.6f  %10.6f  %10.6f" \
                    % (rnum,x1,x2,y1,y2))
            print("           minlevel = %i,  maxlevel = %i" \
                    % (minlevel,maxlevel) \
                    + "  t1 = %s,  t2 = %s" % (f2s(t1),f2s(t2)))
        mapping = {}
        mapping['minlevel'] = minlevel
        mapping['maxlevel'] = maxlevel
        mapping['t1'] = t1
        mapping['t2'] = t2
        mapping['x1'] = x1
        mapping['x2'] = x2
        mapping['y1'] = y1
        mapping['y2'] = y2
        mapping['elev'] = elev
        mapping['name'] = 'Region %i' % rnum
        description = "minlevel = %i, maxlevel = %i\n" % (minlevel,maxlevel) \
            + "  t1 = %s, t2 = %s\n" % (f2s(t1),f2s(t2)) \
            + "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s\n\n" % (f2s(y1),f2s(y2))
        if len(dy_levels) >= minlevel:
            dy = dy_levels[minlevel-1]
            dy_deg,dy_min,dy_sec = deg2dms(dy)
            dy_meters = dy*111e3
            levtext = "Level %s resolution:  \ndy = %g deg, %g min, %g sec \n= %g meters\n" \
                    % (minlevel,dy_deg,dy_min,dy_sec,dy_meters)
            description = description + levtext
        if (maxlevel > minlevel) and (len(dy_levels) >= maxlevel):
            dy = dy_levels[maxlevel-1]
            dy_deg,dy_min,dy_sec = deg2dms(dy)
            dy_meters = dy*111e3
            levtext = "\nLevel %s resolution:  \ndy = %g deg, %g min, %g sec \n= %g meters\n" \
                    % (maxlevel,dy_deg,dy_min,dy_sec,dy_meters)
            description = description + levtext
        mapping['desc'] = description
        mapping['color'] = "FFFFFF"  # white
        mapping['width'] = 3

        region_text = kml_region(mapping)
        kml_text = kml_text + region_text
        if not combined:
            kml_text = kml_text + kml_footer()
            kml_file = open(fname,'w')
            kml_file.write(kml_text)
            kml_file.close()
            if verbose:
                print("Created ",fname)


    try:
        flagregions = rundata.flagregiondata.flagregions
    except:
        flagregions = []  # flagregions not yet in amrclaw

    if len(flagregions)==0 and verbose:
        print("No flagregions found in setrun.py")

    for rnum,flagregion in enumerate(flagregions):

        name = flagregion.name
        #print('+++ flagregion name = ',name)

        if not combined:
            if name == '':
                fname = 'FlagRegion_%s.kml' % str(rnum).zfill(2)
            else:
                fname = name + '.kml'
            kml_text = kml_header(fname)

        #if flagregion.spatial_region is None:
        #    flagregion.read_spatial_region()
        if flagregion.spatial_region_type == 1:
            x1,x2,y1,y2 = flagregion.spatial_region
        else:
            flagregion.read_spatial_region()
            x1,x2,y1,y2 = flagregion.spatial_region.bounding_box()
        minlevel = flagregion.minlevel
        maxlevel = flagregion.maxlevel

        if verbose:
            print("Region ", flagregion.name)

        mapping = {}
        mapping['minlevel'] = flagregion.minlevel
        mapping['maxlevel'] = flagregion.maxlevel
        mapping['t1'] = flagregion.t1
        mapping['t2'] = flagregion.t2
        mapping['x1'] = x1
        mapping['x2'] = x2
        mapping['y1'] = y1
        mapping['y2'] = y2
        mapping['elev'] = elev
        mapping['name'] = flagregion.name
        description = "minlevel = %i, maxlevel = %i\n" \
            % (flagregion.minlevel,flagregion.maxlevel) \
            + "  t1 = %s, t2 = %s\n" % (f2s(flagregion.t1),f2s(flagregion.t2)) \
            + "  Bounding box: \n" \
            + "  x1_bb = %s, x2_bb = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1_bb = %s, y2_bb = %s\n\n" % (f2s(y1),f2s(y2))
        if len(dy_levels) >= minlevel:
            dy = dy_levels[minlevel-1]
            dy_deg,dy_min,dy_sec = deg2dms(dy)
            dy_meters = dy*111e3
            levtext = "Level %s resolution:  \ndy = %g deg, %g min, %g sec \n= %g meters\n" \
                    % (minlevel,dy_deg,dy_min,dy_sec,dy_meters)
            description = description + levtext
        if (maxlevel > minlevel) and (len(dy_levels) >= maxlevel):
            dy = dy_levels[maxlevel-1]
            dy_deg,dy_min,dy_sec = deg2dms(dy)
            dy_meters = dy*111e3
            levtext = "\nLevel %s resolution:  \ndy = %g deg, %g min, %g sec \n= %g meters\n" \
                    % (maxlevel,dy_deg,dy_min,dy_sec,dy_meters)
            description = description + levtext
        mapping['desc'] = description
        mapping['color'] = "00FFFF"  # yellow
        mapping['width'] = 2

        if flagregion.spatial_region_type == 1:
            x1,x2,y1,y2 = flagregion.spatial_region
            x = [x1,x1,x2,x2,x1]
            y = [y1,y2,y2,y1,y1]
        else:
            x,y = flagregion.spatial_region.vertices()

        v = "\n"
        for j in range(len(x)):
            v = v + "%s,%s,%s\n" % (f2s(x[j]),f2s(y[j]),f2s(elev))
        v = v + "%s,%s,%s\n" % (f2s(x[0]),f2s(y[0]),f2s(elev))
        v.replace(' ','')

        region_text = kml_region(mapping, v)

        fname = flagregion.name + '.kml'
        region_text = kml_region(mapping, v)
        kml_text = kml_text + region_text

        if not combined:
            kml_text = kml_text + kml_footer()
            kml_file = open(fname,'w')
            kml_file.write(kml_text)
            kml_file.close()
            if verbose:
                print("Created ",fname)

    if combined:
        fname = fname_combined
        kml_text = kml_text + kml_footer()
        kml_file = open(fname,'w')
        kml_file.write(kml_text)
        kml_file.close()
        if verbose:
            print("Created ",fname)



def line2kml(xy,fname='line.kml',name='line',color='00FFFF',width=3,
             verbose=True):
    """
    Make a KML line with default color yellow.

    :Inputs:

     - *xy* a tuple ((x1,x2),(y1,y2)) (preferred)
            or (x1,x2,y1,y2) (for backward compatibility)
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear on line on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """

    if type(xy[0]) in [tuple,list]:
        x1,x2 = xy[0]
        y1,y2 = xy[1]
    else:
        x1,x2,y1,y2 = xy[0:]

    if verbose:
        print("Line:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s" % (f2s(y1),f2s(y2))
    mapping['color'] = color
    mapping['width'] = width

    region_text = kml_line(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def box2kml(xy,fname=None,name='box',color='FF0000',width=3,verbose=True):
    """
    Make a KML box with default color blue.

    :Inputs:

     - *xy* a tuple ((x1,x2),(y1,y2)) (preferred)
            or (x1,x2,y1,y2) (for backward compatibility)
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear in box on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """

    if fname is None:
        fname = name + '.kml'

    if type(xy[0]) in [tuple,list]:
        x1,x2 = xy[0]
        y1,y2 = xy[1]
    else:
        x1,x2,y1,y2 = xy[0:]

    if verbose:
        print("Box:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %s, x2 = %s\n" % (f2s(x1),f2s(x2)) \
            + "  y1 = %s, y2 = %s" % (f2s(y1),f2s(y2))
    mapping['color'] = color
    mapping['width'] = width

    region_text = kml_region(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def quad2kml(xy,fname=None,name='quad',color='FF0000',width=3,verbose=True):
    """
    Make a KML quadrilateral with default color blue.

    :Inputs:

     - *xy* a tuple ((x1,x2,x3,x4),(y1,y2,y3,y4)) (preferred)
            or (x1,x2,y1,y2,x3,y3,x4,y4) (for backward compatibility)
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear in box on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info

    """

    if fname is None:
        fname = name + '.kml'

    if type(xy[0]) in [tuple,list]:
        x1,x2,x3,x4 = xy[0]
        y1,y2,y3,y4 = xy[1]
    else:
        x1,y1,x2,y2,x3,y3,x4,y4 = xy[0:]

    if verbose:
        print("Quadrilateral:   %10.6f  %10.6f" % (x1,y1))
        print("                 %10.6f  %10.6f" % (x2,y2))
        print("                 %10.6f  %10.6f" % (x3,y3))
        print("                 %10.6f  %10.6f" % (x4,y4))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['x3'] = x3
    mapping['x4'] = x4
    mapping['y1'] = y1
    mapping['y2'] = y2
    mapping['y3'] = y3
    mapping['y4'] = y4
    mapping['elev'] = elev
    mapping['name'] = name
    mapping['desc'] = "  x1 = %s, y1 = %s\n" % (f2s(x1),f2s(y1)) \
            + "  x2 = %s, y2 = %s" % (f2s(x2),f2s(y2)) \
            + "  x3 = %s, y3 = %s" % (f2s(x3),f2s(y3)) \
            + "  x4 = %s, y4 = %s" % (f2s(x4),f2s(y4))
    mapping['color'] = color
    mapping['width'] = 3

    region_text = kml_region(mapping)

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def poly2kml(xy,fname=None,name='poly',color='00FF00', width=3,
             verbose=True, max_vertices_in_description=20):
    """
    Make a KML polygon with default color blue.

    :Inputs:

     - *xy* a tuple (x,y) where x and y are lists of vertices
     - *fname* (str) name of resulting kml file
     - *name* (str) name to appear in box on Google Earth
     - *color* (str) Color in format aabbggrr
     - *width* (str) line width
     - *verbose* (bool) - If *True*, print out info
     - *max_vertices_in_description* (int) - if more than this number
       of vertices, only list number in description box, not all vertices
    """

    if fname is None:
        fname = name + '.kml'

    x,y = xy

    if verbose:
        print("Creating kml for polygon with %i vertices" % len(x))
        if (len(x) <= max_vertices_in_description):
            for j in range(len(x)):
                print("             %10.6f  %10.6f" % (x[j],y[j]))

    elev = 0.
    kml_text = kml_header(fname)

    mapping = {}
    mapping['x'] = x
    mapping['y'] = y
    mapping['elev'] = elev
    mapping['name'] = name
    d = "  Polygon with %i vertices" % len(x)
    if (len(x) <= max_vertices_in_description):
        d = "  x[0] = %s, y[0] = %s\n" % (f2s(x[0]),f2s(y[0]))
        for j in range(1,len(x)):
            d = d + "  x[%i] = %s, y[%i] = %s\n" % (j,f2s(x[j]),j,f2s(y[j]))
    mapping['desc'] = d
    mapping['color'] = color
    mapping['width'] = width

    v = "\n"
    for j in range(len(x)):
        v = v + "%s,%s,%s\n" % (f2s(x[j]),f2s(y[j]),f2s(elev))
    v = v + "%s,%s,%s\n" % (f2s(x[0]),f2s(y[0]),f2s(elev))
    v.replace(' ','')

    region_text = kml_region(mapping, v)
    for j in range(1,len(x)):
        d = d + "  x[%i] = %s, y[%i] = %s" % (j,f2s(x[j]),j,f2s(y[j]))

    kml_text = kml_text + region_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)


def gauges2kml(rundata=None, fname='gauges.kml', verbose=True):

    """

    Create a KML marker for each gauge specified for a GeoClaw run.

    :Inputs:

      - *rundata* - an object of class *ClawRunData* or None

        If *rundata==None*, try to create based on executing function *setrun*
        from the `setrun.py` file in the current directory.

      - *fname* (str) - resulting kml file.

      - *verbose* (bool) - If *True*, print out info about each region found


    :Example:

        >>> from clawpack.geoclaw import kmltools
        >>> kmltools.gauges2kml()

    is equivalent to:

        >>> from clawpack.geoclaw import kmltools
        >>> from setrun import setrun
        >>> rundata = setrun()
        >>> kmltools.gauges2kml(rundata)

    By default this creates a file named *gauges.kml* that can be opened in
    Google Earth.

    """

    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    elev = 0.
    kml_text = kml_header(fname)


    gauges = rundata.gaugedata.gauges
    if len(gauges)==0 and verbose:
        print("No gauges found in setrun.py")


    for rnum,gauge in enumerate(gauges):
        t1,t2 = gauge[3:5]
        x1,y1 = gauge[1:3]
        gaugeno = gauge[0]
        if verbose:
            print("Gauge %i: %s, %s  \n" % (gaugeno,f2s(x1),f2s(y1)) \
                    + "  t1 = %s,  t2 = %s" % (f2s(t1),f2s(t2)))
        mapping = {}
        mapping['gaugeno'] = gaugeno
        mapping['t1'] = t1
        mapping['t2'] = t2
        mapping['x1'] = x1
        mapping['y1'] = y1
        mapping['elev'] = elev
        mapping['name'] = 'Gauge %i' % rnum
        description = "  t1 = %s, t2 = %s\n" % (f2s(t1),f2s(t2)) \
            + "  x1 = %s, y1 = %s\n" % (f2s(x1),f2s(y1))
        mapping['desc'] = description

        gauge_text = kml_gauge(mapping)
        kml_text = kml_text + gauge_text
    kml_text = kml_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()
    if verbose:
        print("Created ",fname)



def kml_header(name='GeoClaw kml file'):
    header = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2"
xmlns:gx="http://www.google.com/kml/ext/2.2">
<Document><name>%s</name>
""" % name
    return header

def kml_footer():
    footer = """
</Document>
</kml>
"""
    return footer


def kml_region(mapping, vertex_text=None):

    if vertex_text is None:
        if 'x3' in mapping:
            # quadrilateral with 4 corners specified
            vertex_text = """
{x1:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y2:.9f},{elev:.9f}
{x3:.9f},{y3:.9f},{elev:.9f}
{x4:.9f},{y4:.9f},{elev:.9f}
{x1:.9f},{y1:.9f},{elev:.9f}
""".format(**mapping).replace(' ','')

        else:
            # rectangle with 2 corners specified
            vertex_text = """
{x1:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y2:.9f},{elev:.9f}
{x1:.9f},{y2:.9f},{elev:.9f}
{x1:.9f},{y1:.9f},{elev:.9f}
""".format(**mapping).replace(' ','')

    mapping['vertices'] = vertex_text
    if len(mapping['color'])==6:
        mapping['color'] = 'FF' + mapping['color']

    kml_text = """
<Style id="Path">
<LineStyle><color>{color:s}</color><width>{width:d}</width></LineStyle>
<PolyStyle><color>00000000</color></PolyStyle>
</Style>
<Placemark><name>{name:s}</name>
<description>{desc:s}</description>
<styleUrl>#Path</styleUrl>
<Polygon>
<tessellate>1</tessellate>
<altitudeMode>clampToGround</altitudeMode>
<outerBoundaryIs><LinearRing><coordinates>
{vertices:s}
</coordinates></LinearRing></outerBoundaryIs>
</Polygon>
</Placemark>
""".format(**mapping)

    return kml_text

def kml_line(mapping):

    if len(mapping['color'])==6:
        mapping['color'] = 'FF' + mapping['color']

    line_text = """
{x1:.9f},{y1:.9f},{elev:.9f}
{x2:.9f},{y2:.9f},{elev:.9f}
""".format(**mapping).replace(' ','')

    mapping['line'] = line_text
    kml_text = """
<Style id="Path">
<LineStyle><color>{color:s}</color><width>{width:d}</width></LineStyle>
<PolyStyle><color>00000000</color></PolyStyle>
</Style>
<Placemark><name>{name:s}</name>
<description>{desc:s}</description>
<styleUrl>#Path</styleUrl>
<LineString>
<tessellate>1</tessellate>
<altitudeMode>clampToGround</altitudeMode>
<coordinates>
{line:s}
</coordinates>
</LineString>
</Placemark>
""".format(**mapping)

    return kml_text

def kml_gauge(mapping):
    gauge_text = "{x1:.9f},{y1:.9f},{elev:.9f}".format(**mapping).replace(' ','')

    mapping['gauge'] = gauge_text

    kml_text = """
<Placemark><name>Gauge {gaugeno:d}</name>
<description>{desc:s}</description>
<styleUrl>#markerstyle</styleUrl>
<Point>
<coordinates>
{gauge:s}
</coordinates>
</Point>
</Placemark>
""".format(**mapping)

    return kml_text



def kml_timespan(t1,t2,event_time=None,tz=None,tscale=1):

    r"""
    Create time strings necessary for sliders in Google Earth.  The time
    span will cover time [t1,t2], with the start of the event given by
    event_time.

    [t1,t2]    : time span,

    event_time : Start of event in UTC :  [Y,M,D,H,M,S], e.g. [2010,2,27,3,34,0]
    tz         : time zone offset to UTC.  e.g. +3 for Chile; -9 for Japan.

    Time span element looks like ::

        <TimeSpan>
          <begin>2010-02-27T06:34:00+03:00</begin>
          <end>2010-02-27T07:04:00+03:00</end>
        </TimeSpan>

    As for how well this handles  Daylight  Savings time, here is what the documentation
    on the Python 'time' module has to say :

    "DST is Daylight Saving Time, an adjustment of the timezone by (usually) one hour
    during part of the year. DST rules are magic (determined by local law) and can
    change from year to year. The C library has a table containing the local rules
    (often it is read from a system file for flexibility) and is the only source of
    True Wisdom in this respect."

    """

    t1 = t1*tscale   # Time converted to seconds
    t2 = t2*tscale

    import time
    # to adjust time from UTC to time in event locale.
    if event_time == None:
        # Use local time.
        starttime = time.mktime(time.localtime())  # seconds UTC
        tz_offset = time.timezone/3600.0   # in seconds
    else:
        ev = tuple(event_time) + (0,0,0)   # Extend to 9 tuple; no DST
        # mktime returns time in seconds + timezone offset, i.e. seconds UTC
        # Subtract out the timezone offset here, since it will get added back
        # in when we do gmtime(starttime + ...) below.
        starttime = time.mktime(ev) - time.timezone
        if tz is None:
            print("===> Time zone offset not defined;  assuming zero offset. " \
                "Set plotdata.kml_tz_offset to define an offset (in hours) from "\
                "UTC (positive west of UTC; negative east of UTC)")
            tz = 0

        tz_offset = tz

    if (tz_offset == None):
        tzstr = "Z"  # no offset; could also just set to "+00:00"
    else:
        # Google Earth will show time slider time in local time, where
        # local + offset = UTC.
        tz_offset = tz_offset*3600.    # Offset in seconds
        tz = time.gmtime(abs(tz_offset))
        if (tz_offset > 0):
            tzstr = time.strftime("+%H:%M",tz)  # Time to UTC
        else:
            tzstr = time.strftime("-%H:%M",tz)

    # Get time strings for start and end of time span
    gbegin = time.gmtime(starttime + t1)
    timestrbegin = "%s%s" % (time.strftime("%Y-%m-%dT%H:%M:%S", gbegin),tzstr)

    gend = time.gmtime(starttime + t2)
    timestrend = "%s%s" % (time.strftime("%Y-%m-%dT%H:%M:%S", gend),tzstr)

    return timestrbegin,timestrend

def topo2kml(topo_file_name, topo_type, color='00FF00'):
    """
    Create a kml file putting a box around the region covered by a topofile.
    Color is green by default.
    """

    import os
    from clawpack.geoclaw import topotools
    topo = topotools.Topography(topo_file_name, topo_type=topo_type)
    topo.read_header()
    xy = topo.extent
    name = os.path.splitext(os.path.split(topo_file_name)[-1])[0]
    file_name = '%s.kml' % name
    box2kml(xy, file_name, name, color)

def dtopo2kml(dtopo_file_name, dtopo_type, color='8888FF'):
    """
    Create a kml file putting a box around the region covered by a dtopofile.
    Color is pink by default.
    """

    import os
    from clawpack.geoclaw import dtopotools
    dtopo = dtopotools.DTopography()
    dtopo.read(dtopo_file_name, dtopo_type)
    x1 = dtopo.x.min()
    x2 = dtopo.x.max()
    y1 = dtopo.y.min()
    y2 = dtopo.y.max()
    xy = (x1,x2,y1,y2)
    name = os.path.splitext(os.path.split(dtopo_file_name)[-1])[0]
    file_name = '%s.kml' % name
    box2kml(xy, file_name, name, color)



def fgmax2kml(rundata=None,fname='fgmax_grids.kml',verbose=True,combined=False):

    """
    Create a KML box for each fgmax grid specified for a GeoClaw run.

    :Inputs:

      - *rundata* - an object of class *ClawRunData* or None

        If *rundata==None*, try to create based on executing function *setrun*
        from the `setrun.py` file in the current directory.

      - *fname* (str) - resulting kml file.

      - *verbose* (bool) - If *True*, print out info about each region found

      - *combined* (bool) - If *True*, combine into single kml file with
        name given by *fname*.  NOT YET IMPLEMENTED.
        If False, *fname* is ignored and individual files are created for
        each fgmax grid.

    """

    from numpy import cos,pi,floor

    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    if combined:
        fname_combined = 'fgmax_grids.kml'
        print('*** combined fgmax kml files not yet supported')
        print('    making a kml file for each fgmax grid')

    fgmax_grids = rundata.fgmax_data.fgmax_grids

    for fg in fgmax_grids:
        fname_root = 'fgmax%s' % str(fg.fgno).zfill(4)
        kml_file = fname_root + '.kml'
        if fg.point_style==1:
            xy = ([fg.x1,fg.x2], [fg.y1,fg.y2])
            line2kml(xy,kml_file, fname_root, color='8888FF', width=2)
        if fg.point_style==2:
            xy = ([fg.x1,fg.x2], [fg.y1,fg.y2])
            box2kml(xy, kml_file, fname_root, color='8888FF')
        elif fg.point_style==3:
            xy = ([fg.x1,fg.x2,fg.x3,fg.x4],
                  [fg.y1,fg.y2,fg.y3,fg.y4])
            poly2kml(xy, kml_file, fname_root, color='8888FF')
        elif fg.point_style==4:
            # points specified by mask in a topo-like file, so plot its extent:
            topo_file_name = fg.xy_fname
            topo_type = 3
            topo2kml(topo_file_name, topo_type, color='8888FF')
        else:
            print('fgmax2kml not yet implemented for point_style = %i' \
                  % fg.point_style)


def fgout2kml(rundata=None,fname='fgout_grids.kml',verbose=True,combined=False):

    """
    Create a KML box for each fgout grid specified for a GeoClaw run.

    :Inputs:

      - *rundata* - an object of class *ClawRunData* or None

        If *rundata==None*, try to create based on executing function *setrun*
        from the `setrun.py` file in the current directory.

      - *fname* (str) - resulting kml file.

      - *verbose* (bool) - If *True*, print out info about each region found

      - *combined* (bool) - If *True*, combine into single kml file with
        name given by *fname*.  NOT YET IMPLEMENTED.
        If False, *fname* is ignored and individual files are created for
        each fgout grid.

    """

    from numpy import cos,pi,floor

    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    if combined:
        fname_combined = 'fgout_grids.kml'
        print('*** combined fgout kml files not yet supported')
        print('    making a kml file for each fgout grid')

    fgout_grids = rundata.fgout_data.fgout_grids

    for fg in fgout_grids:
        fname_root = 'fgout%s' % str(fg.fgno).zfill(4)
        kml_file = fname_root + '.kml'
        xy = ([fg.x1,fg.x2], [fg.y1,fg.y2])
        box2kml(xy, kml_file, fname_root, color='8888FF')


def make_input_data_kmls(rundata=None, combined=False, dtopo_contours=False):
    """
    Produce kml files for the computational domain, all gauges and regions
    specified, and all topo and dtopo files specified in rundata.
    This can be used, e.g. by adding the lines

        from clawpack.geoclaw import kmltools
        kmltools.make_input_data_kmls(rundata)

    to the end of a `setrun.py` file so that `make data` will generate all
    kml files in addition to the `*.data` files.

    Or set *rundata==None*, in which case it will try to generate rundata
    based on executing function *setrun*
    from the `setrun.py` file in the current directory.
    """

    import os
    from clawpack.geoclaw import topotools, dtopotools

    if rundata is None:
        try:
            import setrun
            reload(setrun)
            rundata = setrun.setrun()
        except:
            raise IOError("*** cannot execute setrun file")

    regions2kml(rundata, combined=combined)
    gauges2kml(rundata)
    fgmax2kml(rundata)
    fgout2kml(rundata)

    topofiles = rundata.topo_data.topofiles
    for f in topofiles:
        topo_file_name = f[-1]
        topo_type = f[0]
        topo2kml(topo_file_name, topo_type)

    dtopofiles = rundata.dtopo_data.dtopofiles
    for f in dtopofiles:
        dtopo_file_name = f[-1]
        dtopo_type = f[0]
        dtopo2kml(dtopo_file_name, dtopo_type)
        if dtopo_contours:
            # also make a kmz file showing dtopo contours (takes a bit longer)
            fname_kmz = os.path.split(dtopo_file_name)[-1]
            fname_kmz = os.path.splitext(fname_kmz)[0]  # w/o path or extension
            fname_kmz = fname_kmz + '_contours.kmz'
            dtopo_contours2kmz(dtopo_file_name, dtopo_type=dtopo_type,
                               dZ_interval=1, fname_kmz=fname_kmz,verbose=False)


def pcolorcells_for_kml(X, Y, Z, png_filename=None, dpc=2, max_inches=15.,
                        verbose=True, **kwargs):

    """
    Wraps pcolormesh in a way that a png file is created that can be viewed
    on Google Earth with proper alignment and with sharp grid cell edges.
    Works if X,Y are cell centers or edges, and X,Y can be 2d or 1d arrays.

    X,Y,Z is the data to be plotted.  It is assumed to be finite volume data
    where Z[i,j] is a constant value over a grid cell.

    Internally x,y are defined as 1d arrays since it is assumed the
    grids are Cartesian.

    If the length of the 1d arrays x and y match the dimensions of Z then
    these are assumed to be cell center values. In this case the arrays
    are expanded by one to obtain x_edge, y_edge as edge values,
    as needed for proper alignment.

    If the length of x,y is already one greater than the corresponding
    dimension of Z, then it is assumed that these are already edge values.

    If png_filename is not None then a png file is written with appropriate dpi.

    dpc is the desired "dots per cell", how many pixels to allot to each
    to each grid cell.  This should be an integer to avoid interpolation
    between cells that smears out the cell boundaries in the png file.
    Increasing this will give sharper boundaries but also larger files that
    load more slowly.

    max_inches is the desired size of the longer edge of the figure created.
    This value is not very important unless you want to view the png file
    on a screen outside of Google Earth.  Internally the dimensions of the
    figure `x_inches` and `y_inches` are determined to be consistent
    with the value `dpc` specified and a reasonable value of `dpi` for the
    png file, as described below.

    Internally the value `dpi` (dots per inch) for the png file is
    determined so that it is at least 16 and so that::

        dpi * x_inches = dcp * x_cells
        dpi * y_inches = dcp * y_cells

    where `x_cells`, `y_cells` are the number of cells in each direction.

    `kwargs` are passed to `pcolormesh`, e.g. `cmap` and `norm` are
    generally specified.

    This function returns `fig, ax, png_extent, kml_dpi` so the user can further
    annotate the figure befor saving it as a png file, which should then
    be done with::

        plt.savefig(png_filename, transparent=True, dpi=kml_dpi)

    The `png_extent` is needed in construcing a kml file to display the
    png file on Google Earth, e.g. using the function `png2kml` in this
    module.
    """


    # If X is 2d extract proper 1d slice:
    if X.ndim == 1:
        x = X
    elif X.ndim == 2:
        if X[0,0] == X[0,1]:
            x = X[:,0]
        else:
            x = X[0,:]

    # If Y is 2d extract proper 1d slice:
    if Y.ndim == 1:
        y = Y
    elif Y.ndim == 2:
        if Y[0,0] == Y[0,1]:
            y = Y[:,0]
        else:
            y = Y[0,:]

    dx = x[1]-x[0]
    dy = y[1]-y[0]
    if len(x) == Z.shape[1]:
        # cell centers, so xedge should be expanded by dx/2 on each end:
        xedge = np.arange(x[0]-0.5*dx, x[-1]+dx, dx)
    elif len(x) == Z.shape[1]+1:
        # assume x already contains edge values
        xedge = x
    else:
        raise ValueError('x has unexpected length')

    if len(y) == Z.shape[0]:
        # cell centers, so xedge should be expanded by dx/2 on each end:
        yedge = np.arange(y[0]-0.5*dy, y[-1]+dy, dy)
    elif len(y) == Z.shape[0]+1:
        # assume x already contains edge values
        yedge = y
    else:
        raise ValueError('y has unexpected length')


    x1 = xedge[0];  x2 = xedge[-1]
    y1 = yedge[0];  y2 = yedge[-1]

    # Number of grid cells:
    x_cells = int(round((x2-x1)/dx))
    y_cells = int(round((y2-y1)/dy))
    max_cells = max(x_cells, y_cells)

    dots_per_cell = dpc
    max_dots = dots_per_cell * max_cells

    # determine dots per inch for png file, minimum 64:
    kml_dpi = max(int(round(max_dots / max_inches)), 64)
    dots_x = x_cells * dots_per_cell
    dots_y = y_cells * dots_per_cell

    # determine dimensions for figsize:
    x_inches = dots_x / kml_dpi
    y_inches = dots_y / kml_dpi
    dpc_x = kml_dpi * x_inches / x_cells
    dpc_y = kml_dpi * y_inches / y_cells

    if verbose:
        print('Using kml_dpi = %i,figure size %.6f by %.6f inches' \
                % (kml_dpi,x_inches,y_inches))
        print('Figure has %i by %i grid cells of uniform color' \
                % (x_cells, y_cells))
        print('Dots per cell in x: %.6f,  in y: %.6f' % (dpc_x,dpc_y))
        print('       These should be integers')


    # Create figure of appropriate size and pcolormesh plot
    # with no margins, ticks, or labels:

    fig = plt.figure(figsize=(x_inches,y_inches))
    ax = plt.axes()
    plt.axis('off')
    pc = plt.pcolormesh(xedge, yedge, Z, **kwargs)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    fig.set_size_inches(x_inches,y_inches)
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
                    hspace = 0, wspace = 0)
    plt.margins(0,0)

    if png_filename is not None:
        plt.savefig(png_filename, transparent=True, dpi=kml_dpi)
        if verbose:
            print('Created ',png_filename)

    png_extent = [xedge[0], xedge[-1], yedge[0], yedge[-1]]
    return fig, ax, png_extent, kml_dpi



def kml_png(mapping):
    """
    Create text for a png file overlay
    """

    kml_text = """
<GroundOverlay>
<name>{name:s}</name>
<visibility>1</visibility>
<gx:altitudeMode> clampToSeaFloor </gx:altitudeMode>
<Icon>
  <href>{png_file:s}</href>
</Icon>
<LatLonBox>
  <north>{y2:.9f}</north>
  <south>{y1:.9f}</south>
  <east>{x2:.9f}</east>
  <west>{x1:.9f}</west>
</LatLonBox>
</GroundOverlay>
""".format(**mapping)

    return kml_text

def kml_cb(mapping):
    """
    Create text for a colorbar png file overlay
    """

    kml_text = """
<ScreenOverlay>
  <name>{name:s}</name>
  <Icon>
    <href>{cb_file:s}</href>
  </Icon>
  <overlayXY x="{xfrac:.4f}" xunits="fraction" y="{yfrac:.4f}" yunits="fraction"/>
  <screenXY x="{xfrac:.4f}" xunits="fraction" y="{yfrac:.4f}" yunits="fraction"/>
</ScreenOverlay>
""".format(**mapping)

    return kml_text

radio_style_text = """
<Style id="folderStyle">
<ListStyle>
<listItemType>radioFolder</listItemType>
</ListStyle>
</Style>
<styleUrl>#folderStyle</styleUrl>
"""



def png2kml(extent, png_files, png_names=None, name='png_files', fname=None,
            radio_style=False, cb_files=None, cb_names=None,
            cb_xfracs=None, cb_yfracs=None, verbose=True):
    """
    Create a kml file `fname` linking overlays for each png file in `png_files`.

    `extent` is [x1,x2,y1,y2] specifying where image should be overlaid.

    `png_names`, if present, will give the name for each image for the
    Google Earth menu.

    If `radio_style` is True, set radio buttons so only one can be shown at a
    time, useful for combining plots of different quantities in same file.
    """

    import os

    if fname is None:
        fname = name + '.kml'


    x1,x2,y1,y2 = extent

    if verbose:
        print("Extent:   %10.6f  %10.6f  %10.6f  %10.6f" % (x1,x2,y1,y2))

    kml_text = kml_header(fname) + \
        "<name>%s</name>\n" % name + "<open>1</open>\n"


    mapping = {}
    mapping['x1'] = x1
    mapping['x2'] = x2
    mapping['y1'] = y1
    mapping['y2'] = y2

    for k,png_file in enumerate(png_files):

        mapping['png_file'] = png_file

        try:
            mapping['name'] = png_names[k]
        except:
            mapping['name'] =  'No name'

        kml_text = kml_text + kml_png(mapping)

    if radio_style:
        kml_text = kml_text + radio_style_text

    if cb_files:
        # colorbars
        for k,cb_file in enumerate(cb_files):

            mapping['cb_file'] = cb_file
            try:
                mapping['name'] = cb_names[k]
            except:
                mapping['name'] =  'Colorbar'

            try:
                mapping['xfrac'] = cb_xfracs[k]
            except:
                mapping['xfrac'] = 0.025 + k*0.075
            try:
                mapping['yfrac'] = cb_yfracs[k]
            except:
                mapping['yfrac'] = 0.05

            kml_text = kml_text + kml_cb(mapping)

    kml_text = kml_text + kml_footer()
    kml_file = open(fname,'w')
    kml_file.write(kml_text)
    kml_file.close()

    if verbose:
        print("Created ",fname)


def kml_build_colorbar(cb_filename, cmap, cmin=None, cmax=None,
                       norm=None, label=None, title=None, extend='neither',
                       close_figs=True):

    """
    Make a png file with a colorbar corresponding to cmap, norm.
    cmin, cmax are used only if nrm is not provided.
    """

    import matplotlib as mpl

    fig = plt.figure(figsize=(1.2,3))
    ax1 = fig.add_axes([0.3, 0.075, 0.20, 0.80])
    tick = ax1.yaxis.get_major_ticks()
    plt.tick_params(axis='y', which='major', labelsize=8)

    if norm is None:
        norm = mpl.colors.Normalize(vmin=cmin,vmax=cmax)

    cb1 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,
                                    norm=norm,
                                    extend=extend,
                                    orientation='vertical')

    # make sure ticks appear at lower and upper limits of scale:
    cbticks = cb1.get_ticks()
    cbticks = np.hstack([norm.vmin, cbticks, norm.vmax])
    # remove 2nd and/or next-to-last tick if they are cramped:
    if cbticks[1]-cbticks[0] < 0.25*(cbticks[2]-cbticks[1]):
        cbticks = np.hstack((cbticks[:1], cbticks[2:]))
    if cbticks[-1]-cbticks[-2] < 0.25*(cbticks[-2]-cbticks[-3]):
        cbticks = np.hstack((cbticks[:-2], cbticks[-1:]))
    cb1.set_ticks(cbticks)

    if label:
        cb1.set_label(label)
    if title:
        ax1.set_title(title)


    # This is called from plotpages, in <plotdir>.
    plt.savefig(cb_filename,transparent=False)

    if close_figs:
        plt.close(fig)


def topo2kmz(topo, zlim=(-20,20), mask_outside_zlim=True, sea_level=0.,
             name='topo', force_dry=None, close_figs=True):

    """
    Create kmz file showing onshore and offshore topography as separate layers.
    Currently not very flexible but is useful in quickly creating kmz files
    to open on Google Earth showing the onshore topography and offshore
    bathymetry as two separate layers.  Constant colors on rectangular grid
    cells should be properly rendered to ease exploration of DEMs.
    A colorbar png file is also created and included in the kmz file.

    :Input:
     - *topo* should be a topotools.Topography object,
     - *zlim* is the elevation z limits for choosing the color map
     - *mask_outside_zlim* if True, suppress plotting outsize `zlim`
     - *sea_level* is the break between water and land colors
     - *name* is used in the kml menu and file name
     - *force_dry* is currently not used
     - *close_figs* to close the pyplot figures after making png files

    :Future:
     - If `force_dry` is an array of the same shape as `topo.Z` then another png
       and kml file are created for land that is below `sea_level` but where
       `force_dry = True`.

    """

    import os, glob
    from clawpack.visclaw import colormaps
    import zipfile

    assert force_dry is None, 'force_dry not yet implemented'

    cmap_land = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                         0.25:[0.0,1.0,0.0],
                                          0.5:[0.8,1.0,0.5],
                                          1.0:[0.8,0.5,0.2]})

    cmap_water = colormaps.make_colormap({ 0.0:[0,0,1], 1.:[.8,.8,1]})

    cmap_topo, norm_topo = colormaps.add_colormaps((cmap_land, cmap_water),
                                         data_limits=(zlim[0],zlim[1]),
                                         data_break=sea_level)

    cmap_force_dry = colormaps.make_colormap({ 0.0:[1.0,0.7,0.7], 1.:[1.0,0.7,0.7]})
    cmap_dry, norm_dry = colormaps.add_colormaps((cmap_land, cmap_force_dry),
                                         data_limits=(zlim[0],zlim[1]),
                                         data_break=sea_level)

    if mask_outside_zlim:
        Z = ma.masked_where(np.logical_or(topo.Z<zlim[0], topo.Z>zlim[1]), topo.Z)
        cbar_extend = 'neither'
    else:
        Z = topo.Z
        cbar_extend = 'both'

    kml_dir = 'kmlfiles_%s' % name
    os.system('mkdir -p %s' % kml_dir)
    print('Will put png and kml files in %s' % kml_dir)

    Z_land = ma.masked_where(Z<sea_level, Z)
    png_filename = '%s/%s_land.png' % (kml_dir, name)
    fig,ax,png_extent,kml_dpi = pcolorcells_for_kml(topo.X, topo.Y,
                                         Z_land, png_filename=png_filename,
                                         dpc=2, cmap=cmap_topo, norm=norm_topo)
    if close_figs:
        plt.close(fig)

    Z_water = ma.masked_where(Z>=sea_level, Z)
    png_filename = '%s/%s_water.png' % (kml_dir, name)
    fig,ax,png_extent,kml_dpi = pcolorcells_for_kml(topo.X, topo.Y,
                                         Z_water, png_filename=png_filename,
                                         dpc=2, cmap=cmap_topo, norm=norm_topo)

    if close_figs:
        plt.close(fig)

    kml_build_colorbar('%s/colorbar.png' % kml_dir, cmap_topo,
                                cmin=zlim[0], cmax=zlim[1], label='meters',
                                title='topo', extend=cbar_extend,
                                close_figs=close_figs)


    png_files=['%s_water.png' % name, '%s_land.png' % name]
    png_names=['%s_water' % name,'%s_land' % name]
    cb_files = ['colorbar.png']
    cb_names = ['colorbar_topo']

    name = '%s_topo' % name
    fname = os.path.join(kml_dir, name+'.kml')
    png2kml(png_extent, png_files=png_files, png_names=png_names,
                     name=name, fname=fname,
                     radio_style=False,
                     cb_files=cb_files, cb_names=cb_names)

    savedir = os.getcwd()
    os.chdir(kml_dir)
    files = glob.glob('*.kml') + glob.glob('*.png')
    print('kmz file will include:')
    for file in files:
        print('    %s' % os.path.split(file)[-1])

    fname_kmz = 'topo_%s.kmz' % name
    with zipfile.ZipFile(fname_kmz, 'w') as zip:
        for file in files:
            zip.write(file)
        print('Created %s' % os.path.abspath(fname_kmz))
    os.chdir(savedir)

def transect2kml(xtrans, ytrans, fname='transect.kml'):
    """
    Create a kml file for points with long,lat specified by xtrans,ytrans.
    Label each point with number and coordinates.
    Adjust longitudes in coordinates so they are all between -180 and 180
    to display properly in Google Earth (but labels are original values).
    """

    with open(fname,'w') as kml_file:

        kml_file.write("""<?xml version="1.0" encoding="UTF-8"?>
        <kml xmlns="http://www.opengis.net/kml/2.2"
        xmlns:gx="http://www.google.com/kml/ext/2.2">
        <Document><name>transect_points</name>

        <Style id="Red">
        <IconStyle><color>FF0000FF</color></IconStyle>
        </Style>

        <Style id="Yellow">
        <IconStyle><color>FF00FFFF</color></IconStyle>
        </Style>
        """)

        for k in range(len(xtrans)):
            x = xtrans[k]
            y = ytrans[k]
            if x < -180:
                xge = x+360
            elif x > 180:
                xge = x-360
            else:
                xge = x
            style_id = 'Red'
            name = 'Transect point %i' % k
            kml_file.write("""
            <Placemark><name>%s</name>
            <description>%.5f, %.5f</description>
            <styleUrl>#%s</styleUrl>
            <Point>
            <coordinates>
            %.9f, %.9f, 0.0
            </coordinates>
            </Point>
            </Placemark>
            """ % (name,x,y,style_id,xge,y))

        kml_file.write("\n</Document>\n</kml>")

    print('Created ', fname)

def dtopo_contours2kmz(dtopofiles, dtopo_type=3, dZ_interval=1, dZmax=40,
                       text_label=True, text_x=None, text_y=None,
                       fname_kmz=None, verbose=True):

    """
    Create dtopo_contours.kmz file containing contour plots of the dtopo
    deformations with radio buttons to select which one to show.
    (Or dtopofiles can be a string for a single dtopofile.)

    To show multiple dtopofiles, they must all have the same extent.

    :Inputs:

     - *dtopofiles* (str or list of str): single or list of dtopofile paths
     - *dZ_interval* (float) the interval in meters between contours.
     - *dZmax* (float) max for contour levels shown (and -dZmax is the min)
     - *text_label* Text label to add to plots,
        If text_label is a string this will be added as a text label,
        If text_label == True a standard label will be added with the
        dopofile name, dZ_interval, and the max, min values of dZ.
     - *text_x, text_y* are the location of the text,
        or if None then the mean of dtopofile.X and dtopofile.Y are used.

    :Future:
     - Could also add pcolor plot of deformation to kmz file.
    """
    from clawpack.geoclaw import dtopotools

    import os,sys, zipfile, shutil

    kml_dir = 'temp_dtopos_kml'
    os.system('mkdir -p %s' % kml_dir)
    events = []
    png_files = []

    if type(dtopofiles) == str:
        dtopofiles = [dtopofiles]  # if only one passed in

    previous_png_extent = None

    for dtopofile in dtopofiles:
        dtopofile = os.path.abspath(dtopofile)
        event = os.path.splitext(os.path.split(dtopofile)[-1])[0]
        events.append(event)
        if verbose:
            print('Making contours for event = ',event)
            print('  contour interval = %gm' % dZ_interval)

        dtopo = dtopotools.DTopography(dtopofile, dtopo_type=dtopo_type)

        # first make empty plot of right dimensions:
        Zm = ma.masked_where(dtopo.Y < 90, dtopo.Y)  # all masked
        fig,ax,png_extent,kml_dpi = pcolorcells_for_kml(dtopo.X, dtopo.Y, Zm,
                                                 png_filename=None, dpc=4,
                                                 verbose=verbose)

        if previous_png_extent is not None:
            assert png_extent == previous_png_extent, \
                    '*** dtopofiles have different extents\n' \
                  + '*** cannot make a single kmz file displaying them'
        previous_png_extent = png_extent

        # now add contour plot:
        clines = np.arange(dZ_interval,dZmax,dZ_interval)
        lw = 2.
        ax.contour(dtopo.X, dtopo.Y, dtopo.dZ[-1,:,:], clines, colors='r',
                   linestyles='-', linewidths=lw)
        clines = np.arange(-dZmax,0,dZ_interval)
        ax.contour(dtopo.X, dtopo.Y, dtopo.dZ[-1,:,:], clines, colors='c',
                   linestyles='-', linewidths=lw)

        if text_label:
            if text_x is None:
                text_x = dtopo.x.mean()
            if text_y is None:
                text_y = dtopo.y.mean()

            dZ_min = dtopo.dZ[-1,:,:].min()
            dZ_max = dtopo.dZ[-1,:,:].max()
            if type(text_label) == str:
                plt.text(text_x, text_y, text_label, fontsize=15,color='yellow')
            else:
                plt.text(text_x, text_y,'%s\ndZ_min = %.1fm, dZ_max = %.1fm' \
                        % (event,dZ_min,dZ_max) \
                        +  '\ndZ_interval = %.1f m' % dZ_interval, \
                     fontsize=15,color='yellow')

        png_filename = '%s_contours.png' % event
        png_files.append(png_filename)
        png_filename = os.path.join(kml_dir, png_filename)
        plt.savefig(png_filename, transparent=True, dpi=kml_dpi)
        plt.close(fig)

    if fname_kmz is None:
        fname_kmz = 'dtopo_contours.kmz'
    path_kmz = os.path.abspath('./%s' % fname_kmz)  # path in this directory

    # move to kml_dir for making kml files:
    savedir = os.getcwd()
    os.chdir(kml_dir)
    png_names = events
    fname = 'dtopo_contours.kml'
    name = 'dtopo_contours'
    png2kml(png_extent, png_files=png_files, png_names=png_names,
                     radio_style=True, name=name, fname=fname, verbose=verbose)

    files = [fname] + ['%s_contours.png' % event for event in events]
    if 0:
        print('kmz file will include:')
        for file in files:
            print('    ', file)

    with zipfile.ZipFile(path_kmz, 'w') as zip:
        for file in files:
            zip.write(file)
        print('Created %s' % fname_kmz)
    os.chdir(savedir)

    if verbose:
        print('Removing ',kml_dir)
    shutil.rmtree(kml_dir)
