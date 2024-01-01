import numpy as np
from datetime import datetime
import os
from dataclasses import dataclass


@dataclass
# Creates a dataclass for holding the OWI data
class StormData:
    iLat: int  # number of latitude points
    iLong: int  # number of longitude points
    DX: float  # resolution in x direction
    DY: float  # resolution in y direction
    SWLat: float  # Initial Latitude point in SW corner
    SWLon: float  # Initial Longitude point in SW corner
    DT: datetime  # datestamp of current wind/pressure array
    matrix: list  # placeholder for wind or pressure array

    def __post_init__(self):
        # Put everything in correct format
        self.iLat = int(self.iLat)
        self.iLong = int(self.iLong)
        self.DX = float(self.DX)
        self.DY = float(self.DY)
        self.SWLat = float(self.SWLat)
        self.SWLon = float(self.SWLon)
        self.DT = datetime.strptime(self.DT, '%Y%m%d%H%M')

def read_oceanweather(path, file_ext):
    import data_storms
    import re
    """
        Reads in Oceanweather files and puts them into a dataclass for ease of data cleanup
        """
    subset = None
    all_data = []
    fh = path + f'.{file_ext}'
    # Open file and use regex matching for parsing data
    with open(fh, 'rt') as f:
        input = f.readlines()
        print(type(input), len(input))
        for line in input:
            if not line.startswith('Oceanweather'):  # Skip the file header
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
    t0 = None
    time_array = []
    total_time = data[-1].DT - data[0].DT
    num_steps = int(round((total_time / len(data)).total_seconds() / 60))
    for idx, d in enumerate(data):
        if not t0:
            t0 = d.DT
        t = d.DT
        seconds_from_start = (t - t0).total_seconds()
        time_array.append(seconds_from_start)
    return time_array


def get_coordinate_arrays(data):
    lat = [data.SWLat + i * data.DY for i in range(data.iLat)]
    lon = [data.SWLon + i * data.DX for i in range(data.iLong)]
    return lat, lon

def arrange_data(data):
    data_list = [item for sublist in data.matrix for item in sublist]
    return data_list


def process_data(data, start_idx=0):
    values = arrange_data(data)
    d = np.empty(shape=(data.iLong, data.iLat))
    for j in range(data.iLat):
        for i in range(data.iLong):
            d[i,j] = values[start_idx + j * data.iLong + i]
    return d.T # Transpose to fit into the correct format for geoclaw

def find_eye(pressure_data):
    # Calculate the location of the eye for all pressure arrays
    # Return None if the min pressure > 10000 because that
    # means the storm's eye is outside the domain

    # if pressure_data is in mbars or Pa
    if len(str(int(pressure_data[0].max()))) == 4:
        min_p = 990.00
    else:
        min_p = 99000.00
    eye_locs = [np.unravel_index(np.argmin(pressure), pressure.shape)
                if np.min(pressure) < min_p else None for pressure in pressure_data]


    return eye_locs


def find_last_closed_isobar(lat, lon, pressure_array, eye):
    import matplotlib.pyplot as plt

    # Use matplotlib contour function to find the closed contours of the pressure array
    lev = [950, 960, 970, 980, 990, 1000, 1005, 1010]
    contour = plt.contour(lon, lat, pressure_array, levels=lev)

    target_value = 1000
    tolerance = 1.0

    #Extract contour indices
    contour_inds = contour.collections
    target_contour_index = -1
    target_contour_indices = []

    for idx, contour_set in enumerate(contour_inds):
        for contour_line in contour_set.get_paths():
            x_values = contour_line.vertices[:,0]
            y_values = contour_line.vertices[:,1]
            if abs(contour.levels[idx] - target_value) < tolerance:
                target_contour_index = idx
                target_contour_indices = list(zip(x_values, y_values))
    # # Variables to store information about the largest closed contour
    # max_closed_contour_index = -1
    # max_closed_contour_area = 0
    # max_closed_inds = []
    # max_distance = 350000 # KM
    #
    # for idx, contour_set in enumerate(contour_inds):
    #     for contour_line in contour_set.get_paths():
    #         x_values = contour_line.vertices[:,0]
    #         y_values = contour_line.vertices[:,1]
    #
    #         # Check if the contour is closed (first_point == last_point)
    #         if x_values[0] == x_values[-1] and y_values[0] == y_values[-1]:
    #             # Calcualte teh area of the closed contour
    #             contour_area = np.abs(np.trapz(y_values, x_values))
    #
    #             # Check if the current closed contour is larger than the last
    #             if contour_area > max_closed_contour_area:
    #                 max_closed_contour_area = contour_area
    #                 max_closed_contour_index = idx
    #                 max_closed_inds = list(zip(x_values, y_values))
    #
    #
    # if max_closed_contour_index == -1:
    #     for idx, contour_set in enumerate(contour_inds):
    #         for contour_line in contour_set.get_paths():
    #             x_values = contour_line.vertices[:,0]
    #             y_values = contour_line.vertices[:,1]
    #
    #             distance_to_point = [haversine(lat, lon, eye[0], eye[1])
    #                                  for lat, lon in zip(y_values, x_values)]
    #             print(distance_to_point)
    #             if np.min(distance_to_point) < max_distance:
    #                 contour_area = np.abs(np.trapz(y_values, x_values))
    #
    #                 if contour_area > max_closed_contour_area:
    #                     max_closed_contour_area = contour_area
    #                     max_closed_contour_index = idx
    #                     max_closed_inds = list(zip(x_values, y_values))

    return np.array(target_contour_indices)



def haversine(lat1, lon1, lat2, lon2):
    import math
    # Radius of the Earth in kilometers
    R = 6371000.0

    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(abs(lon1))
    lat2 = math.radians(lat2)
    lon2 = math.radians(abs(lon2))

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c

    return distance

def pressure_radius(lat, lon, pressure, eye):
    # Find the indices in all directions of the first time the pressure
    # changes to atmospheric, calculate the average to return the storm radius
    dist = []
    if eye:
        lat1 = lat[eye[0]]
        lon1 = lon[eye[1]]
        indices = find_last_closed_isobar(lat, lon, pressure, (lat1, lon1))
        for pair in indices:
            lat2 = pair[1]
            lon2 = pair[0]
            dist.append(haversine(lat1, lon1, lat2, lon2))
        return np.average(dist)
    else:
        return None


def calc_rmw(lat, lon, data, eye_loc):
    if eye_loc:
        y, x = np.unravel_index(data.argmax(), data.shape)
        lat1 = lat[eye_loc[0]]
        lon1 = lon[eye_loc[1]]
        lat2 = lat[y]
        lon2 = lon[x]
        dist = haversine(lat1, lon1, lat2, lon2)
        return dist
    else:
        return None





