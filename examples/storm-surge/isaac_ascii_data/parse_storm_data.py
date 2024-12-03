from clawpack.geoclaw.surge import storm
import os

"""
This script will convert a wind and pressure field pair of files into a netcdf
for use in geoclaw's data derived storm objects

The ability to run both netcdf files and the original .win and .pre files is now implemented as well
In the setrun.py set the storm_specification_type as -2 for netcdf or -3 for ascii data derived inputs
"""

storm_path = os.environ['CLAW']
storm_file = os.path.join(storm_path, 'isaac_data')
# Instantiate the class
storm = storm.DataDerivedStorms(storm_file, wind_file_ext='WIN', pressure_file_ext='PRE')
storm.parse_data(landfall_time='2012-08-29T00:00:00')

output_name = os.path.join(storm_path,'isaac_data_derived')
storm.write_data_derived(output_name)
