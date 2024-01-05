import numpy as np
from datetime import datetime
from dataclasses import dataclass

"""
A suite of helper functions for parsing data derived storms in the Oceanweather fixed width format
and converting to a DataDerivedStorm class for GeoClaw
"""
@dataclass
# Creates a dataclass for holding the OWI data
class StormData:
    i_lat: int  # number of latitude points
    i_long: int  # number of longitude points
    dx: float  # resolution in x direction
    dy: float  # resolution in y direction
    sw_lat: float  # Initial Latitude point in SW corner
    sw_lon: float  # Initial Longitude point in SW corner
    dt: datetime  # datestamp of current wind/pressure array
    matrix: list  # placeholder for wind or pressure array

    def __post_init__(self):
        # Put everything in correct format
        self.i_lat = int(self.i_lat) # Number of latitude points
        self.i_long = int(self.i_long) # Number of longitude points
        self.dx = float(self.dx) # resolution in x direction
        self.dy = float(self.dy) # resolution in y direction
        self.sw_lat = float(self.sw_lat) # Initial latitude point in sw corner
        self.sw_lon = float(self.sw_lon) # Initial longitude point in sw corner

        # Check if DT is not already a datetime object before conversion
        if not isinstance(self.dt, datetime):
            self.dt = datetime.strptime(self.dt, '%Y%m%d%H%M')

def read_oceanweather(path, file_ext):
    import re
    """
        Reads in Oceanweather files and puts them into a dataclass for ease of data cleanup
        
        :param path: The path to the file.
        :param file_ext: The file extension.
        :return: A list of StormData objects.
        """
    subset = None
    all_data = []

    # Open file and use regex matching for parsing data
    with open("{path}.{file_ext}", 'rt') as f:
        input = f.readlines()
        print(type(input), len(input))
        for line in input:
            if not line.startswith('Oceanweather'):  # Skip the Oceanweather file header
                # Find the header lines containing this pattern of characters
                # Example from file: iLat= 105iLong=  97DX=0.2500DY=0.2500SWLat=22.00000SWLon=-82.0000DT=200007121200
                # \w+: any unicode string
                # \s*: any repetition of whitespace characters
                # \d+: decimal digit of length +=1
                # \.?: matches anything but a new line in minimal fashion
                # \d*: decimal digit with +=0 repetitions
                header = re.findall("\\w+=-?\\s*\\d+\\.?\\d*", line)
                if header:
                    if subset:
                        # put data into dataclass
                        storm_data = StormData(**subset)
                        all_data.append(storm_data)

                    # Split apart the header data into separate values rather than the string
                    subset = {
                        x.replace(' ', '').split('=')[0]: x.replace(' ', '').split('=')[1]
                        for x in header
                    }
                    subset["matrix"] = []
                else:
                    nums = list(map(float, line.split()))
                    subset["matrix"].append(nums)
                storm_data = StormData(**subset)
        all_data.append(storm_data)

    return all_data

def time_steps(data):
    """
    Calculate the timesteps for geoclaw in seconds given the start and end
    times in the header of the wind and pre files. Returns an array of time steps
    with 0 being at the start of the data
    :param data: wind or pressure data read from read_oceanweather
    :return: array of time steps with total length = length of data
    """
    start_time = None
    time_array = []
    for idx, d in enumerate(data):
        if not start_time:
            start_time = d.dt
        t = d.dt
        seconds_from_start = (t - start_time).total_seconds()
        time_array.append(seconds_from_start)
    return time_array


def get_coordinate_arrays(data):
    """
    Creates a list of latitude and longitude points given the starting location in the sw corner
    and uses the resolution and number of points per array
    :param data: StormData
    :return: list
    """
    lat = [data.sw_lat + i * data.dy for i in range(data.i_lat)]
    lon = [data.sw_lon + i * data.dx for i in range(data.i_long)]
    return lat, lon

def arrange_data(data):
    """
    Iterates over the entire matrix for a single timestep and formats the data
    into a single array rather than a list of lists
    :param data: StormData
    :return: list
    """
    data_list = [item for sublist in data.matrix for item in sublist]
    return data_list


def process_data(data, start_idx=0):
    """
    Process wind and pressure data into a 2d array

    :param data: StormData
    :param start_idx: starting point depending on u or v direction
    :return: Array
    """
    # Flatten list of lists into a single list
    values = arrange_data(data)

    #Fill 2d array with values
    d = np.empty(shape=(data.i_long, data.i_lat))
    for j in range(data.i_lat):
        for i in range(data.i_long):
            d[i,j] = values[start_idx + j * data.i_long + i]
    return d.T # Transpose to fit into the correct format for geoclaw






