#!/usr/bin/env python

"""

Classes representing parameters for 1D GeoClaw runs

:Classes:

 - GeoClawData1D
 - GaugeData1D

:Constants:

 - Rearth - Radius of earth in meters
 - DEG2RAD factor to convert degrees to radians
 - RAD2DEG factor to convert radians to degrees

"""

import os
import numpy
import clawpack.clawutil.data

# Radius of earth in meters.
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
LAT2METER = Rearth * DEG2RAD

class GeoClawData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D Geoclaw data object

    """
    def __init__(self):
        super(GeoClawData1D,self).__init__()

        # GeoClaw physics parameters
        self.add_attribute('gravity',9.8)
        self.add_attribute('earth_radius',Rearth)
        self.add_attribute('coordinate_system',1)
        self.add_attribute('friction_forcing',True)
        self.add_attribute('friction_coefficient',0.025)

        # GeoClaw algorithm parameters
        self.add_attribute('dry_tolerance',1e-3)
        self.add_attribute('friction_depth',1.0e6)
        self.add_attribute('sea_level',0.0)


    def write(self,data_source='setrun.py'):

        self.open_data_file('geoclaw.data',data_source)

        self.data_write('gravity')
        self.data_write('earth_radius')
        self.data_write('coordinate_system')
        self.data_write('sea_level')

        friction = self.friction_forcing
        if isinstance(self.friction_forcing,bool):
            if self.friction_forcing:
                friction = 1
            else:
                friction = 0
        elif isinstance(self.friction_forcing,str):
            if self.friction_forcing in ['Manning','manning','MANNING']:
                friction = 1
            elif self.friction_forcing in ['Coulomb','coulomb','COULOMB']:
                friction = 2
            else:
                friction = 0
        self.friction_forcing = friction


        self.data_write('friction_forcing')
        self.data_write('friction_coefficient')
        self.data_write('friction_depth')
        self.data_write('dry_tolerance')

        self.close_data_file()

#  Gauge data object

#  Gauge data object
class GaugeData1D(clawpack.clawutil.data.ClawData):
    r"""
     Gauge data object for 1d.
     input specs for gauges are in 1d in setrun.py...output is like that of
      2d amr (with level=1 and y=0) so that same reading/plotting tools can be used.
    """

    @property
    def gauge_numbers(self):
        if len(self.gauges) == 1:
            return [self.gauges[0][0]]
        else:
            return [gauge[0] for gauge in self.gauges]

    def __init__(self, num_dim=2):
        super(GaugeData1D,self).__init__()

        self.add_attribute('num_dim',num_dim)
        self.add_attribute('gauges',[])

    def __str__(self):
        output = "Gauges: %s\n" % len(self.gauges)
        for gauge in self.gauges:
            output = "\t".join((output,"%4i:" % gauge[0]))
            output = " ".join((output,"%19.10e" % gauge[1]))
            output = " ".join((output,"%17.10e" % gauge[2]))
            output = " ".join((output,"%13.6e\n" % gauge[3]))
        return output

    def write(self,out_file='gauges.data',data_source='setrun.py'):
        r"""Write out gague information data file."""


        # Check to make sure we have only unique gauge numebrs
        if len(self.gauges) > 0:
            if len(self.gauge_numbers) != len(set(self.gauge_numbers)):
                raise Exception("Non unique gauge numbers specified.")

        # Write out gauge data file
        self.open_data_file(out_file,data_source)
        self.data_write(name='ngauges',value=len(self.gauges))
        for gauge in self.gauges:
            self._out_file.write("%4i %19.10e %19.10e %13.6e  %13.6e\n" % tuple(gauge))
        self.close_data_file()

    def read(self,data_path="./",file_name='gauges.data'):
        r"""Read gauge data file"""
        path = os.path.join(data_path, file_name)
        gauge_file = open(path,'r')

        # Read past comments and blank lines
        header_lines = 0
        ignore_lines = True
        while ignore_lines:
            line = gauge_file.readline()
            if line[0] == "#" or len(line.strip()) == 0:
                header_lines += 1
            else:
                break

        # Read number of gauges, should be line that was last read in
        num_gauges = int(line.split()[0])

        # Read in each gauge line
        for n in xrange(num_gauges):
            line = gauge_file.readline().split()
            self.gauges.append([int(line[0]),float(line[1]),float(line[2]),
                                                              float(line[3]),float(line[4])])

        gauge_file.close()


