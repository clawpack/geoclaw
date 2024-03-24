from clawpack.geoclaw.surge import storm
import os


storm_path = os.environ['CLAW']
storm_file = os.path.join(storm_path, 'isaac_data')
# Instantiate the class
storm = storm.DataDerivedStorms(storm_file, wind_file_ext='WIN', pressure_file_ext='PRE')
storm.parse_data()

output_name = os.path.join(storm_path,'isaac_data_derived')
storm.write_data_derived(output_name)