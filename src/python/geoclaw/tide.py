#!/usr/bin/env python

"""
GeoClaw util Module  `$CLAW/geoclaw/src/python/geoclaw/tide.py`

Module provides provides tide prediction functions.

:Functions:

 - retrieve_constituents - retrieves harmonic constituents from NOAA gauge station
 - retrieve_water_levels - retrieves observed water levels from NOAA's API
 - retrieve_predicted_tide - retrieves predicted tide from NOAA's API
 - datum_value - retrieves datum value for desired datum reference
 - predict_tide - predicts tide with Pytides
 - datetimes - prepares a collection of datetimes from beginning to end dates
 - detide - detides observed water levels
 - surge - predicts surge at NOAA gauge station
"""

from __future__ import absolute_import
from __future__ import print_function
from collections.abc import Iterable
from collections import OrderedDict, namedtuple
from scipy.optimize import leastsq, fsolve
from itertools import takewhile, count
from datetime import datetime, timedelta
from functools import reduce
from six.moves.urllib.parse import urlencode
from six.moves.urllib.request import urlopen

try:
    from itertools import izip, ifilter
except ImportError: #Python3
    izip = zip
    ifilter = filter
    
try:
    import requests
    import json
    import string
    import lxml.html as lh
    import pandas as pd
    import operator as op
    import numpy as np
    import io
    import os
    import os.path
except ImportError as e:
    print(e)

#%env CLAW=/Users/jonathansocoy/clawpack

d2r, r2d = np.pi/180.0, 180.0/np.pi

NOAA_API_URL = 'https://tidesandcurrents.noaa.gov/api/datagetter'
NOAA_home = 'https://tidesandcurrents.noaa.gov/harcon.html'

######################## Tide Prediction Functions ########################

def retrieve_constituents(station, time_zone='GMT', units='meters', cache_dir=None,
                         verbose=True):
                         
    """Fetch constituent data for given NOAA tide station.
    By default, retrieved data is cached in the geoclaw scratch directory
    located at:
        $CLAW/geoclaw/scratch
    :Required Arguments:
      - station (string): 7 character station ID
    :Optional Arguments:
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
      - cache_dir (string): alternative directory to use for caching data
      - verbose (bool): whether to output informational messages
    :Returns:
      - constituent_dict (dictionary): dictionary of tidal constituents for NOAA gauge station
    """
                        
    def units_num(units):
        if (units == 'meters'):
            return 0
        elif (time_zone == 'feet'):
            return 1
    
    def time_zone_num(time_zone):
        if (time_zone == 'GMT'):
            return 0
        elif (time_zone == 'Local'):
            return 1
        
    def get_noaa_params(station, time_zone, units):
        noaa_params = {
            'unit': units_num(units),
            'timezone': time_zone_num(time_zone),
            'id': station
        }
        return noaa_params
    
    # use geoclaw scratch directory for caching by default
    if cache_dir is None:
        if 'CLAW' not in os.environ:
            raise ValueError('CLAW environment variable not set')
        claw_dir = os.environ['CLAW']
        cache_dir = os.path.join(claw_dir, 'geoclaw', 'scratch')    #### cache_dir
 
    def get_cache_path(station, time_zone, units):
        filename = '{}_{}_constituents'.format(time_zone, units)
        abs_cache_dir = os.path.abspath(cache_dir)
        return os.path.join(abs_cache_dir, 'constituents', station, filename)
    
    def save_to_cache(cache_path, data):
        # make parent directories if they do not exist
        parent_dir = os.path.dirname(cache_path)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)
            
        component_array = pd.DataFrame(data)
        component_array.to_csv(cache_path, index=False)
            
    def parse(cache_path):
        # read data into structured array, skipping header row if present
        data = pd.read_csv(cache_path)
        component_array = pd.DataFrame(data)
        component_dict = component_array.to_dict(orient='list')
        return component_dict

    #Requests URL
    def fetch_data(station, time_zone, units):
        noaa_params = get_noaa_params(station, time_zone, units)
        cache_path = get_cache_path(station, time_zone, units)
        
        # use cached data if available
        if os.path.exists(cache_path):
            if verbose:
                print('Using cached constituent data for station {}'.format(station))
            return parse(cache_path)
            

        # otherwise, retrieve data from NOAA and cache it
        if verbose:
            print('Fetching constituent data from NOAA for station {}'.format(station))
  
        #Forms URL
        url = '{}?{}'.format(NOAA_home, urlencode(noaa_params))

        page = requests.get(url)
        doc = lh.fromstring(page.content)
        tr_elements = doc.xpath('//tr')
        col = [((t.text_content(),[])) for t in tr_elements[0]]
        for j in range(1, len(tr_elements)):
            T, i = tr_elements[j], 0
            for t in T.iterchildren():
                col[i][1].append(t.text_content())
                i+=1
                
        constituent_dict = {title:column for (title,column) in col}

        # if there were no errors, then cache response
        save_to_cache(cache_path, constituent_dict)
                
        return constituent_dict


    try:
        constituents = fetch_data(station, time_zone, units)
         
    except:
        print('*** Fetching NOAA Constituents failed, returning None')
        constituents = None
    
    return constituents


def fetch_noaa_tide_data(station, begin_date, end_date, datum='MTL', time_zone='GMT', units='metric', cache_dir=None, verbose=True):
    """Fetch water levels and tide predictions at given NOAA tide station.
    The data is returned in 6 minute intervals between the specified begin and
    end dates/times.  A complete specification of the NOAA CO-OPS API for Data
    Retrieval used to fetch the data can be found at:
        https://tidesandcurrents.noaa.gov/api/
    By default, retrieved data is cached in the geoclaw scratch directory
    located at:
        $CLAW/geoclaw/scratch
    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of retrieval
      - end_date (datetime): end of date/time range of retrieval
    :Optional Arguments:
      - datum (string): see NOAA API documentation for possible values
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
      - cache_dir (string): alternative directory to use for caching data
      - verbose (bool): whether to output informational messages
    :Returns:
      - date_time (numpy.ndarray): times corresponding to retrieved data
      - water_level (numpy.ndarray): preliminary or verified water levels
      - prediction (numpy.ndarray): tide predictions
    """
    # use geoclaw scratch directory for caching by default
    if cache_dir is None:
        if 'CLAW' not in os.environ:
            raise ValueError('CLAW environment variable not set')
        claw_dir = os.environ['CLAW']
        cache_dir = os.path.join(claw_dir, 'geoclaw', 'scratch')

    def fetch(product, expected_header, col_idx, col_types):
        noaa_params = get_noaa_params(product)
        cache_path = get_cache_path(product)

        # use cached data if available
        if os.path.exists(cache_path):
            if verbose:
                print('Using cached {} data for station {}'.format(
                    product, station))
            return parse(cache_path, col_idx, col_types, header=True)

        # otherwise, retrieve data from NOAA and cache it
        if verbose:
            print('Fetching {} data from NOAA for station {}'.format(
                product, station))
        full_url = '{}?{}'.format(NOAA_API_URL, urlencode(noaa_params))
        with urlopen(full_url) as response:
            text = response.read().decode('utf-8')
            with io.StringIO(text) as data:
                # ensure that received header is correct
                header = data.readline().strip()
                if header != expected_header or 'Error' in text:
                    # if not, response contains error message
                    raise ValueError(text)

                # if there were no errors, then cache response
                save_to_cache(cache_path, text)

                return parse(data, col_idx, col_types, header=False)

    def get_noaa_params(product):
        noaa_date_fmt = '%Y%m%d %H:%M'
        noaa_params = {
            'product': product,
            'application': 'Clawpack',
            'format': 'csv',
            'station': station,
            'begin_date': begin_date.strftime(noaa_date_fmt),
            'end_date': end_date.strftime(noaa_date_fmt),
            'time_zone': time_zone,
            'datum': datum,
            'units': units
        }
        return noaa_params

    def get_cache_path(product):
        cache_date_fmt = '%Y%m%d%H%M'
        dates = '{}_{}'.format(begin_date.strftime(cache_date_fmt),
                               end_date.strftime(cache_date_fmt))
        filename = '{}_{}_{}'.format(time_zone, datum, units)
        abs_cache_dir = os.path.abspath(cache_dir)
        return os.path.join(abs_cache_dir, product, station, dates, filename)

    def save_to_cache(cache_path, data):
        # make parent directories if they do not exist
        parent_dir = os.path.dirname(cache_path)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

        # write data to cache file
        with open(cache_path, 'w') as cache_file:
            cache_file.write(data)

    def parse(data, col_idx, col_types, header):
        # read data into structured array, skipping header row if present
        a = np.genfromtxt(data, usecols=col_idx, dtype=col_types,
                             skip_header=int(header), delimiter=',',
                             missing_values='')

        # return tuple of columns
        return tuple(a[col] for col in a.dtype.names)

    # only need first two columns of data; first column contains date/time,
    # and second column contains corresponding value
    col_idx = (0, 1)
    col_types = 'datetime64[m], float'

    # fetch water levels and tide predictions
    try:
        date_time, water_level = fetch(
            'water_level', 'Date Time, Water Level, Sigma, O or I (for verified), F, R, L, Quality',
            col_idx, col_types)
    except:
        print('*** Fetching water_level failed, returning None')
        date_time = None
        water_level = None

    try:
        date_time2, prediction = fetch('predictions', 'Date Time, Prediction',
                                       col_idx, col_types)
        if date_time is None:
            date_time = date_time2

    except:
        print('*** Fetching prediction failed, returning None')
        date_time2 = None
        prediction = None

    # ensure that date/time ranges are the same
    if (date_time is not None) and (date_time2 is not None):
        if not np.array_equal(date_time, date_time2):
            raise ValueError('Received data for different times')

    return date_time, water_level, prediction
    
    
def datum_value(station, datum, time_zone='GMT', units='metric'):
    """Fetch datum value for given NOAA tide station.
    :Required Arguments:
      - station (string): 7 character station ID
      - datum (string): MSL, MTL
    :Optional Arguments:
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
    :Returns:
      - datum_value (float): value for requested datum reference
    """
    def get_noaa_params(station, time_zone, units):
        noaa_params = {
            'product': 'datums',
            'units': units,
            'time_zone': time_zone,
            'station': station,
            'application': 'Clawpack',
            'format':'json'
        }
        return noaa_params
    
    #Scrapes MTL/MSL Datum Value
    def get_datum(station, time_zone, units):
        noaa_params = get_noaa_params(station, time_zone, units)
        url = '{}?{}'.format(NOAA_API_URL, urlencode(noaa_params))
        page_data = requests.get(url)
        data = page_data.json()['datums']
        datum_value = [d['v'] for d in data if d['n'] == datum]
        return datum_value
    
    try:
        datum_value = float(get_datum(station, time_zone, units)[0])
    except:
        print('*** Fetching datum value failed, returning None')
        datum_value = None
        
    return datum_value

     
def predict_tide(station, begin_date, end_date, datum='MTL', time_zone='GMT', units='meters'):
    """Fetch datum value for given NOAA tide station.
    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of prediction
      - end_date (datetime): end of date/time range of prediction
    :Optional Arguments:
      - datum (string): MTL for tide prediction
      - time_zone (string): see NOAA API documentation for possible values
      - units (string): see NOAA API documentation for possible values
    :Returns:
      - heights (float): tide heights
    """
    #These are the NOAA constituents, in the order presented on NOAA's website.
    constituents = [c for c in noaa if c != _Z0]
    noaa_values = retrieve_constituents(station, time_zone, units)
    noaa_amplitudes = [float(amplitude) for amplitude in noaa_values['Amplitude']]
    noaa_phases = [float(phases) for phases in noaa_values['Phase']] 
    
    #We can add a constant offset - set to MTL
#    MTL = datum_value(args[0], 'MTL')
    desired_datum = datum_value(station, datum)
    MSL = datum_value(station, 'MSL')
    offset = MSL - desired_datum
    constituents.append(_Z0)
    noaa_phases.append(0)
    noaa_amplitudes.append(offset)
       
    #Build the model
    assert(len(constituents) == len(noaa_phases) == len(noaa_amplitudes))
    model = np.zeros(len(constituents), dtype = Tide.dtype)
    model['constituent'] = constituents
    model['amplitude'] = noaa_amplitudes
    model['phase'] = noaa_phases
    tide = Tide(model = model, radians = False)
    
    #Time Calculations
    delta = (end_date-begin_date)/timedelta(hours=1) + .1
    times = Tide._times(begin_date, np.arange(0, delta, .1))
    
    #Height Calculations
    heights_arrays = [tide.at([times[i]]) for i in range(len(times))]
    heights = [val for sublist in heights_arrays for val in sublist]
 
    return heights


def datetimes(begin_date, end_date):
    #Time Calculations
    delta = (end_date-begin_date)/timedelta(hours=1) + .1
    times = Tide._times(begin_date, np.arange(0, delta, .1))
    return times


def detide(NOAA_observed_water_level, predicted_tide):
    # NOAA observed water level - predicted tide 
    return [(NOAA_observed_water_level[i] - predicted_tide[i]) for i in range(len(NOAA_observed_water_level))]

#Surge Implementation
def surge(station, beg_date, end_date, landfall_date):
    """Fetch datum value for given NOAA tide station.
    :Required Arguments:
      - station (string): 7 character station ID
      - begin_date (datetime): start of date/time range of prediction
      - end_date (datetime): end of date/time range of prediction
      - landfall_date (datetime): approximate time of landfall for reference
    :Optional Arguments:
      - datum (string): MTL for tide prediction and retrieval
      - time_zone (string): see NOAA API documentation for possible values
    :Returns:
      - times (float): times with landfall event as reference
      - surge (float): surge heights
    """
    predicted_tide = predict_tide(station, beg_date, end_date)
    NOAA_times, NOAA_observed_water_level, NOAA_predicted_tide = fetch_noaa_tide_data(station, beg_date, end_date)

    #detides NOAA observed water levels with predicted tide
    surge = detide(NOAA_observed_water_level, predicted_tide)
    #modifies NOAA times to datetimes
    times = [((pd.to_datetime(time).to_pydatetime())-landfall_date)/timedelta(days=1) for time in NOAA_times]
    
    return times, surge



######################## Nodal Corrections ########################

def f_unity(a):
	return 1.0

#Schureman equations 73, 65
def f_Mm(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = (2/3.0 - np.sin(omega)**2)*(1 - 3/2.0 * np.sin(i)**2)
	return (2/3.0 - np.sin(I)**2) / mean

#Schureman equations 74, 66
def f_Mf(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(omega)**2 * np.cos(0.5*i)**4
	return np.sin(I)**2 / mean

#Schureman equations 75, 67
def f_O1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(omega) * np.cos(0.5*omega)**2 * np.cos(0.5*i)**4
	return (np.sin(I) * np.cos(0.5*I)**2) / mean

#Schureman equations 76, 68
def f_J1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(2*omega) * (1-3/2.0 * np.sin(i)**2)
	return np.sin(2*I) / mean

#Schureman equations 77, 69
def f_OO1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.sin(omega) * np.sin(0.5*omega)**2 * np.cos(0.5*i)**4
	return np.sin(I) * np.sin(0.5*I)**2 / mean

#Schureman equations 78, 70
def f_M2(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	mean = np.cos(0.5*omega)**4 * np.cos(0.5*i)**4
	return np.cos(0.5*I)**4 / mean

#Schureman equations 227, 226, 68
#Should probably eventually include the derivations of the magic numbers (0.5023 etc).
def f_K1(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	nu = d2r*a['nu'].value
	sin2Icosnu_mean = np.sin(2*omega) * (1-3/2.0 * np.sin(i)**2)
	mean = 0.5023*sin2Icosnu_mean + 0.1681
	return (0.2523*np.sin(2*I)**2 + 0.1689*np.sin(2*I)*np.cos(nu)+0.0283)**(0.5) / mean

#Schureman equations 215, 213, 204
#It can be (and has been) confirmed that the exponent for R_a reads 1/2 via Schureman Table 7
def f_L2(a):
	P = d2r*a['P'].value
	I = d2r*a['I'].value
	R_a_inv = (1 - 12*np.tan(0.5*I)**2 * np.cos(2*P)+36*np.tan(0.5*I)**4)**(0.5)
	return f_M2(a) * R_a_inv

#Schureman equations 235, 234, 71
#Again, magic numbers
def f_K2(a):
	omega = d2r*a['omega'].value
	i = d2r*a['i'].value
	I = d2r*a['I'].value
	nu = d2r*a['nu'].value
	sinsqIcos2nu_mean = np.sin(omega)**2 * (1-3/2.0 * np.sin(i)**2)
	mean = 0.5023*sinsqIcos2nu_mean + 0.0365
	return (0.2533*np.sin(I)**4 + 0.0367*np.sin(I)**2 *np.cos(2*nu)+0.0013)**(0.5) / mean

#Schureman equations 206, 207, 195
def f_M1(a):
	P = d2r*a['P'].value
	I = d2r*a['I'].value
	Q_a_inv = (0.25 + 1.5*np.cos(I)*np.cos(2*P)*np.cos(0.5*I)**(-0.5) + 2.25*np.cos(I)**2 * np.cos(0.5*I)**(-4))**(0.5)
	return f_O1(a) * Q_a_inv

#See e.g. Schureman equation 149
def f_Modd(a, n):
	return f_M2(a) ** (n / 2.0)

#Node factors u, see Table 2 of Schureman.

def u_zero(a):
	return 0.0

def u_Mf(a):
	return -2.0 * a['xi'].value

def u_O1(a):
	return 2.0 * a['xi'].value - a['nu'].value

def u_J1(a):
	return -a['nu'].value

def u_OO1(a):
	return -2.0 * a['xi'].value - a['nu'].value

def u_M2(a):
	return 2.0 * a['xi'].value - 2.0 * a['nu'].value

def u_K1(a):
	return -a['nup'].value

#Schureman 214
def u_L2(a):
	I = d2r*a['I'].value
	P = d2r*a['P'].value
	R = r2d*np.arctan(np.sin(2*P)/(1/6.0 * np.tan(0.5*I) **(-2) -np.cos(2*P)))
	return 2.0 * a['xi'].value - 2.0 * a['nu'].value - R

def u_K2(a):
	return -2.0 * a['nupp'].value

#Schureman 202
def u_M1(a):
	I = d2r*a['I'].value
	P = d2r*a['P'].value
	Q = r2d*np.arctan((5*np.cos(I)-1)/(7*np.cos(I)+1)*np.tan(P))
	return a['xi'].value - a['nu'].value + Q

def u_Modd(a, n):
	return n/2.0 * u_M2(a)



######################## Constituents ########################


class BaseConstituent(object):
	xdo_int = {
		'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8, 'I': 9,
		'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 'O': 15, 'P': 16, 'Q': 17,
		'R': -8, 'S': -7, 'T': -6, 'U': -5, 'V': -4, 'W': -3, 'X': -2, 'Y': -1,
		'Z': 0
	}

	int_xdo = {v:k for k, v in xdo_int.items()}

	def __init__(self, name, xdo='', coefficients=[], u=u_zero, f=f_unity):
		if xdo == '':
			self.coefficients = np.array(coefficients)
		else:
			self.coefficients = np.array(self.xdo_to_coefficients(xdo))
		self.name = name
		self.u = u
		self.f = f

	def xdo_to_coefficients(self, xdo):
		return [self.xdo_int[l.upper()] for l in xdo if l in string.ascii_letters]

	def coefficients_to_xdo(self, coefficients):
		return ''.join([self.int_xdo[c] for c in coefficients])

	def V(self, astro):
		return np.dot(self.coefficients, self.astro_values(astro))

	def xdo(self):
		return self.coefficients_to_xdo(self.coefficients)

	def speed(self, a):
		return np.dot(self.coefficients, self.astro_speeds(a))

	def astro_xdo(self, a):
		return [a['T+h-s'], a['s'], a['h'], a['p'], a['N'], a['pp'], a['90']]

	def astro_speeds(self, a):
		return np.array([each.speed for each in self.astro_xdo(a)])

	def astro_values(self, a):
		return np.array([each.value for each in self.astro_xdo(a)])

	#Consider two out of phase constituents which travel at the same speed to
	#be identical
	def __eq__(self, c):
		return np.all(self.coefficients[:-1] == c.coefficients[:-1])

	def __hash__(self):
		return hash(tuple(self.coefficients[:-1]))

class CompoundConstituent(BaseConstituent):

	def __init__(self, members = [], **kwargs):
		self.members = members

		if 'u' not in kwargs:
			kwargs['u'] = self.u
		if 'f' not in kwargs:
			kwargs['f'] = self.f

		super(CompoundConstituent,self).__init__(**kwargs)

		self.coefficients = reduce(op.add,[c.coefficients * n for (c,n) in members])

	def speed(self, a):
		return reduce(op.add, [n * c.speed(a) for (c,n) in self.members])

	def V(self, a):
		return reduce(op.add, [n * c.V(a) for (c,n) in self.members])

	def u(self, a):
		return reduce(op.add, [n * c.u(a) for (c,n) in self.members])

	def f(self, a):
		return reduce(op.mul, [c.f(a) ** abs(n) for (c,n) in self.members])

###### Base Constituents
#Long Term
_Z0      = BaseConstituent(name = 'Z0',      xdo = 'Z ZZZ ZZZ', u = u_zero, f = f_unity)
_Sa      = BaseConstituent(name = 'Sa',      xdo = 'Z ZAZ ZZZ', u = u_zero, f = f_unity)
_Ssa     = BaseConstituent(name = 'Ssa',     xdo = 'Z ZBZ ZZZ', u = u_zero, f = f_unity)
_Mm      = BaseConstituent(name = 'Mm',      xdo = 'Z AZY ZZZ', u = u_zero, f = f_Mm)
_Mf      = BaseConstituent(name = 'Mf',      xdo = 'Z BZZ ZZZ', u = u_Mf, f = f_Mf)

#Diurnals
_Q1      = BaseConstituent(name = 'Q1',      xdo = 'A XZA ZZA', u = u_O1, f = f_O1)
_O1      = BaseConstituent(name = 'O1',      xdo = 'A YZZ ZZA', u = u_O1, f = f_O1)
_K1      = BaseConstituent(name = 'K1',      xdo = 'A AZZ ZZY', u = u_K1, f = f_K1)
_J1      = BaseConstituent(name = 'J1',      xdo = 'A BZY ZZY', u = u_J1, f = f_J1)

#M1 is a tricky business for reasons of convention, rather than theory.  The
#reasons for this are best summarised by Schureman paragraphs 126, 127 and in
#the comments found in congen_input.txt of xtides, so I won't go over all this
#again here.

_M1      = BaseConstituent(name = 'M1',      xdo = 'A ZZZ ZZA', u = u_M1, f = f_M1)
_P1      = BaseConstituent(name = 'P1',      xdo = 'A AXZ ZZA', u = u_zero, f = f_unity)
_S1      = BaseConstituent(name = 'S1',      xdo = 'A AYZ ZZZ', u = u_zero, f = f_unity)
_OO1     = BaseConstituent(name = 'OO1',     xdo = 'A CZZ ZZY', u = u_OO1, f = f_OO1)

#Semi-Diurnals
_2N2     = BaseConstituent(name = '2N2',     xdo = 'B XZB ZZZ', u = u_M2, f = f_M2)
_N2      = BaseConstituent(name = 'N2',      xdo = 'B YZA ZZZ', u = u_M2, f = f_M2)
_nu2     = BaseConstituent(name = 'nu2',     xdo = 'B YBY ZZZ', u = u_M2, f = f_M2)
_M2      = BaseConstituent(name = 'M2',      xdo = 'B ZZZ ZZZ', u = u_M2, f = f_M2)
_lambda2 = BaseConstituent(name = 'lambda2', xdo = 'B AXA ZZB', u = u_M2, f = f_M2)
_L2      = BaseConstituent(name = 'L2',      xdo = 'B AZY ZZB', u = u_L2, f = f_L2)
_T2      = BaseConstituent(name = 'T2',      xdo = 'B BWZ ZAZ', u = u_zero, f = f_unity)
_S2      = BaseConstituent(name = 'S2',      xdo = 'B BXZ ZZZ', u = u_zero, f = f_unity)
_R2      = BaseConstituent(name = 'R2',      xdo = 'B BYZ ZYB', u = u_zero, f = f_unity)
_K2      = BaseConstituent(name = 'K2',      xdo = 'B BZZ ZZZ', u = u_K2, f = f_K2)

#Third-Diurnals
_M3      = BaseConstituent(name = 'M3',      xdo = 'C ZZZ ZZZ', u = lambda a: u_Modd(a,3), f = lambda a: f_Modd(a,3))

###### Compound Constituents
#Long Term
_MSF     = CompoundConstituent(name = 'MSF',  members = [(_S2, 1), (_M2, -1)])

#Diurnal
_2Q1     = CompoundConstituent(name = '2Q1',  members = [(_N2, 1), (_J1, -1)])
_rho1    = CompoundConstituent(name = 'rho1', members = [(_nu2, 1), (_K1, -1)])

#Semi-Diurnal

_mu2     = CompoundConstituent(name = 'mu2',  members = [(_M2, 2), (_S2, -1)]) #2MS2
_2SM2    = CompoundConstituent(name = '2SM2', members = [(_S2, 2), (_M2, -1)])

#Third-Diurnal
_2MK3    = CompoundConstituent(name = '2MK3', members = [(_M2, 1), (_O1, 1)])
_MK3     = CompoundConstituent(name = 'MK3',  members = [(_M2, 1), (_K1, 1)])

#Quarter-Diurnal
_MN4     = CompoundConstituent(name = 'MN4',  members = [(_M2, 1), (_N2, 1)])
_M4      = CompoundConstituent(name = 'M4',   members = [(_M2, 2)])
_MS4     = CompoundConstituent(name = 'MS4',  members = [(_M2, 1), (_S2, 1)])
_S4      = CompoundConstituent(name = 'S4',   members = [(_S2, 2)])

#Sixth-Diurnal
_M6      = CompoundConstituent(name = 'M6',   members = [(_M2, 3)])
_S6      = CompoundConstituent(name = 'S6',   members = [(_S2, 3)])

#Eighth-Diurnals
_M8      = CompoundConstituent(name = 'M8',   members = [(_M2, 4)])


noaa = [
	_M2, _S2, _N2, _K1, _M4, _O1, _M6, _MK3, _S4, _MN4, _nu2, _S6, _mu2,
	_2N2, _OO1, _lambda2, _S1, _M1, _J1, _Mm, _Ssa, _Sa, _MSF, _Mf,
	_rho1, _Q1, _T2, _R2, _2Q1, _P1, _2SM2, _M3, _L2, _2MK3, _K2,
	_M8, _MS4
]


####################### Tide ########################


class Tide(object):
	dtype = np.dtype([
		('constituent', object),
		('amplitude', float),
		('phase', float)])

	def __init__(
			self,
			constituents = None,
			amplitudes = None,
			phases = None,
			model = None,
			radians = False
		):
		"""
		Initialise a tidal model. Provide constituents, amplitudes and phases OR a model.
		Arguments:
		constituents -- list of constituents used in the model.
		amplitudes -- list of amplitudes corresponding to constituents
		phases -- list of phases corresponding to constituents
		model -- an ndarray of type Tide.dtype representing the constituents, amplitudes and phases.
		radians -- boolean representing whether phases are in radians (default False)
		"""
		if None not in [constituents, amplitudes, phases]:
			if len(constituents) == len(amplitudes) == len(phases):
				model = np.zeros(len(phases), dtype=Tide.dtype)
				model['constituent'] = np.array(constituents)
				model['amplitude'] = np.array(amplitudes)
				model['phase'] = np.array(phases)
			else:
				raise ValueError("Constituents, amplitudes and phases should all be arrays of equal length.")
		elif model is not None:
			if not model.dtype == Tide.dtype:
				raise ValueError("Model must be a numpy array with dtype == Tide.dtype")
		else:
			raise ValueError("Must be initialised with constituents, amplitudes and phases; or a model.")
		if radians:
			model['phase'] = r2d*model['phase']
		self.model = model[:]
		self.normalize()

	def prepare(self, *args, **kwargs):
		return Tide._prepare(self.model['constituent'], *args, **kwargs)

	@staticmethod
	def _prepare(constituents, t0, t = None, radians = True):
		"""
		Return constituent speed and equilibrium argument at a given time, and constituent node factors at given times.
		Arguments:
		constituents -- list of constituents to prepare
		t0 -- time at which to evaluate speed and equilibrium argument for each constituent
		t -- list of times at which to evaluate node factors for each constituent (default: t0)
		radians -- whether to return the angular arguments in radians or degrees (default: True)
		"""
		#The equilibrium argument is constant and taken at the beginning of the
		#time series (t0).  The speed of the equilibrium argument changes very
		#slowly, so again we take it to be constant over any length of data. The
		#node factors change more rapidly.
		if isinstance(t0, Iterable):
			t0 = t0[0]
		if t is None:
			t = [t0]
		if not isinstance(t, Iterable):
			t = [t]
		a0 = astro(t0)
		a = [astro(t_i) for t_i in t]

		#For convenience give u, V0 (but not speed!) in [0, 360)
		V0 = np.array([c.V(a0) for c in constituents])[:, np.newaxis]
		speed = np.array([c.speed(a0) for c in constituents])[:, np.newaxis]
		u = [np.mod(np.array([c.u(a_i) for c in constituents])[:, np.newaxis], 360.0)
			 for a_i in a]
		f = [np.mod(np.array([c.f(a_i) for c in constituents])[:, np.newaxis], 360.0)
			 for a_i in a]

		if radians:
			speed = d2r*speed
			V0 = d2r*V0
			u = [d2r*each for each in u]
		return speed, u, f, V0

	def at(self, t):
		"""
		Return the modelled tidal height at given times.
		Arguments:
		t -- array of times at which to evaluate the tidal height
		"""
		t0 = t[0]
		hours = self._hours(t0, t)
		partition = 240.0
		t = self._partition(hours, partition)
		times = self._times(t0, [(i + 0.5)*partition for i in range(len(t))])
		speed, u, f, V0 = self.prepare(t0, times, radians = True)
		H = self.model['amplitude'][:, np.newaxis]
		p = d2r*self.model['phase'][:, np.newaxis]

		return np.concatenate([
			Tide._tidal_series(t_i, H, p, speed, u_i, f_i, V0)
			for t_i, u_i, f_i in izip(t, u, f)
		])

	def highs(self, *args):
		"""
		Generator yielding only the high tides.
		Arguments:
		see Tide.extrema()
		"""
		for t in ifilter(lambda e: e[2] == 'H', self.extrema(*args)):
			yield t

	def lows(self, *args):
		"""
		Generator yielding only the low tides.
		Arguments:
		see Tide.extrema()
		"""
		for t in ifilter(lambda e: e[2] == 'L', self.extrema(*args)):
			yield t

	def form_number(self):
		"""
		Returns the model's form number, a helpful heuristic for classifying tides.
		"""
		k1, o1, m2, s2 = (
			np.extract(self.model['constituent'] == c, self.model['amplitude'])
			for c in [_K1, _O1, _M2, _S2]
		)
		return (k1+o1)/(m2+s2)

	def classify(self):
		"""
		Classify the tide according to its form number
		"""
		form = self.form_number()
		if 0 <= form <= 0.25:
			return 'semidiurnal'
		elif 0.25 < form <= 1.5:
			return 'mixed (semidiurnal)'
		elif 1.5 < form <= 3.0:
			return 'mixed (diurnal)'
		else:
			return 'diurnal'

	def extrema(self, t0, t1 = None, partition = 2400.0):
		"""
		A generator for high and low tides.
		Arguments:
		t0 -- time after which extrema are sought
		t1 -- optional time before which extrema are sought (if not given, the generator is infinite)
		partition -- number of hours for which we consider the node factors to be constant (default: 2400.0)
		"""
		if t1:
			#yield from in python 3.4
			for e in takewhile(lambda t: t[0] < t1, self.extrema(t0)):
				yield e
		else:
			#We assume that extrema are separated by at least delta hours
			delta = np.amin([
				90.0 / c.speed(astro(t0)) for c in self.model['constituent']
				if not c.speed(astro(t0)) == 0
			])
			#We search for stationary points from offset hours before t0 to
			#ensure we find any which might occur very soon after t0.
			offset = 24.0
			partitions = (
				Tide._times(t0, i*partition) for i in count()), (Tide._times(t0, i*partition) for i in count(1)
			)

			#We'll overestimate to be on the safe side;
			#values outside (start,end) won't get yielded.
			interval_count = int(np.ceil((partition + offset) / delta)) + 1
			amplitude = self.model['amplitude'][:, np.newaxis]
			phase     = d2r*self.model['phase'][:, np.newaxis]

			for start, end in izip(*partitions):
				speed, [u], [f], V0 = self.prepare(start, Tide._times(start, 0.5*partition))
				#These derivatives don't include the time dependence of u or f,
				#but these change slowly.
				def d(t):
					return np.sum(-speed*amplitude*f*np.sin(speed*t + (V0 + u) - phase), axis=0)
				def d2(t):
					return np.sum(-speed**2.0 * amplitude*f*np.cos(speed*t + (V0 + u) - phase), axis=0)
				#We'll overestimate to be on the safe side;
				#values outside (start,end) won't get yielded.
				intervals = (
					delta * i -offset for i in range(interval_count)), (delta*(i+1) - offset for i in range(interval_count)
				)
				for a, b in izip(*intervals):
					if d(a)*d(b) < 0:
						extrema = fsolve(d, (a + b) / 2.0, fprime = d2)[0]
						time = Tide._times(start, extrema)
						[height] = self.at([time])
						hilo = 'H' if d2(extrema) < 0 else 'L'
						if start < time < end:
							yield (time, height, hilo)

	@staticmethod
	def _hours(t0, t):
		"""
		Return the hourly offset(s) of a (list of) time from a given time.
		Arguments:
		t0 -- time from which offsets are sought
		t -- times to find hourly offsets from t0.
		"""
		if not isinstance(t, Iterable):
			return Tide._hours(t0, [t])[0]
		elif isinstance(t[0], datetime):
			return np.array([(ti-t0).total_seconds() / 3600.0 for ti in t])
		else:
			return t

	@staticmethod
	def _partition(hours, partition = 3600.0):
		"""
		Partition a sorted list of numbers (or in this case hours).
		Arguments:
		hours -- sorted ndarray of hours.
		partition -- maximum partition length (default: 3600.0)
		"""
		partition = float(partition)
		relative = hours - hours[0]
		total_partitions = np.ceil(relative[-1] / partition + 10*np.finfo(np.float).eps).astype('int')
		return [hours[np.floor(np.divide(relative, partition)) == i] for i in range(total_partitions)]

	@staticmethod
	def _times(t0, hours):
		"""
		Return a (list of) datetime(s) given an initial time and an (list of) hourly offset(s).
		Arguments:
		t0 -- initial time
		hours -- hourly offsets from t0
		"""
		if not isinstance(hours, Iterable):
			return Tide._times(t0, [hours])[0]
		elif not isinstance(hours[0], datetime):
			return np.array([t0 + timedelta(hours=h) for h in hours])
		else:
			return np.array(hours)

	@staticmethod
	def _tidal_series(t, amplitude, phase, speed, u, f, V0):
		return np.sum(amplitude*f*np.cos(speed*t + (V0 + u) - phase), axis=0)

	def normalize(self):
		"""
		Adapt self.model so that amplitudes are positive and phases are in [0,360) as per convention
		"""
		for i, (_, amplitude, phase) in enumerate(self.model):
			if amplitude < 0:
				self.model['amplitude'][i] = -amplitude
				self.model['phase'][i] = phase + 180.0
			self.model['phase'][i] = np.mod(self.model['phase'][i], 360.0)

	@classmethod
	def decompose(
			cls,
			heights,
			t            = None,
			t0           = None,
			interval     = None,
			constituents = noaa,
			initial      = None,
			n_period     = 2,
			callback     = None,
			full_output  = False
		):
		"""
		Return an instance of Tide which has been fitted to a series of tidal observations.
		Arguments:
		It is not necessary to provide t0 or interval if t is provided.
		heights -- ndarray of tidal observation heights
		t -- ndarray of tidal observation times
		t0 -- datetime representing the time at which heights[0] was recorded
		interval -- hourly interval between readings
		constituents -- list of constituents to use in the fit (default: noaa)
		initial -- optional Tide instance to use as first guess for least squares solver
		n_period -- only include constituents which complete at least this many periods (default: 2)
		callback -- optional function to be called at each iteration of the solver
		full_output -- whether to return the output of scipy's leastsq solver (default: False)
		"""
		if t is not None:
			if isinstance(t[0], datetime):
				hours = Tide._hours(t[0], t)
				t0 = t[0]
			elif t0 is not None:
				hours = t
			else:
				raise ValueError("t can be an array of datetimes, or an array "
				                 "of hours since t0 in which case t0 must be "
				                 "specified.")
		elif None not in [t0, interval]:
			hours = np.arange(len(heights)) * interval
		else:
			raise ValueError("Must provide t(datetimes), or t(hours) and "
			                 "t0(datetime), or interval(hours) and t0(datetime) "
			                 "so that each height can be identified with an "
			                 "instant in time.")

		#Remove duplicate constituents (those which travel at exactly the same
		#speed, irrespective of phase)
		constituents = list(OrderedDict.fromkeys(constituents))

		#No need for least squares to find the mean water level constituent z0,
		#work relative to mean
		constituents = [c for c in constituents if not c == _Z0]
		z0 = np.mean(heights)
		heights = heights - z0

		#Only analyse frequencies which complete at least n_period cycles over
		#the data period.
		constituents = [
			c for c in constituents
			if 360.0 * n_period < hours[-1] * c.speed(astro(t0))
		]
		n = len(constituents)

		sort = np.argsort(hours)
		hours = hours[sort]
		heights = heights[sort]

		#We partition our time/height data into intervals over which we consider
		#the values of u and f to assume a constant value (that is, their true
		#value at the midpoint of the interval).  Constituent
		#speeds change much more slowly than the node factors, so we will
		#consider these constant and equal to their speed at t0, regardless of
		#the length of the time series.

		partition = 240.0

		t     = Tide._partition(hours, partition)
		times = Tide._times(t0, [(i + 0.5)*partition for i in range(len(t))])

		speed, u, f, V0 = Tide._prepare(constituents, t0, times, radians = True)

		#Residual to be minimised by variation of parameters (amplitudes, phases)
		def residual(hp):
			H, p = hp[:n, np.newaxis], hp[n:, np.newaxis]
			s = np.concatenate([
				Tide._tidal_series(t_i, H, p, speed, u_i, f_i, V0)
				for t_i, u_i, f_i in izip(t, u, f)
			])
			res = heights - s
			if callback:
				callback(res)
			return res

		#Analytic Jacobian of the residual - this makes solving significantly
		#faster than just using gradient approximation, especially with many
		#measurements / constituents.
		def D_residual(hp):
			H, p = hp[:n, np.newaxis], hp[n:, np.newaxis]
			ds_dH = np.concatenate([
				f_i*np.cos(speed*t_i+u_i+V0-p)
				for t_i, u_i, f_i in izip(t, u, f)],
				axis = 1)

			ds_dp = np.concatenate([
				H*f_i*np.sin(speed*t_i+u_i+V0-p)
				for t_i, u_i, f_i in izip(t, u, f)],
				axis = 1)

			return np.append(-ds_dH, -ds_dp, axis=0)

		#Initial guess for solver, haven't done any analysis on this since the
		#solver seems to converge well regardless of the initial guess We do
		#however scale the initial amplitude guess with some measure of the
		#variation
		amplitudes = np.ones(n) * (np.sqrt(np.dot(heights, heights)) / len(heights))
		phases     = np.ones(n)

		if initial:
			for (c0, amplitude, phase) in initial.model:
				for i, c in enumerate(constituents):
					if c0 == c:
						amplitudes[i] = amplitude
						phases[i] = d2r*phase

		initial = np.append(amplitudes, phases)

		lsq = leastsq(residual, initial, Dfun=D_residual, col_deriv=True, ftol=1e-7)

		model = np.zeros(1+n, dtype=cls.dtype)
		model[0] = (_Z0, z0, 0)
		model[1:]['constituent'] = constituents[:]
		model[1:]['amplitude'] = lsq[0][:n]
		model[1:]['phase'] = lsq[0][n:]

		if full_output:
			return cls(model = model, radians = True), lsq
		return cls(model = model, radians = True)
    

################# Astronomical Values #######################


#Convert a sexagesimal angle into decimal degrees
def s2d(degrees, arcmins = 0, arcsecs = 0, mas = 0, muas = 0):
	return (
			degrees
			+ (arcmins /  60.0)
			+ (arcsecs / (60.0*60.0))
			+ (mas	   / (60.0*60.0*1e3))
			+ (muas    / (60.0*60.0*1e6))
	)

#Evaluate a polynomial at argument
def polynomial(coefficients, argument):
	return sum([c * (argument ** i) for i,c in enumerate(coefficients)])

#Evaluate the first derivative of a polynomial at argument
def d_polynomial(coefficients, argument):
	return sum([c * i * (argument ** (i-1)) for i,c in enumerate(coefficients)])

#Meeus formula 11.1
def T(t):
	return (JD(t) - 2451545.0)/36525

#Meeus formula 7.1
def JD(t):
	Y, M = t.year, t.month
	D = (
		t.day
		+ t.hour / (24.0)
		+ t.minute / (24.0*60.0)
		+ t.second / (24.0*60.0*60.0)
		+ t.microsecond / (24.0 * 60.0 * 60.0 * 1e6)
	)
	if M <= 2:
		Y = Y - 1
		M = M + 12
	A = np.floor(Y / 100.0)
	B = 2 - A + np.floor(A / 4.0)
	return np.floor(365.25*(Y+4716)) + np.floor(30.6001*(M+1)) + D + B - 1524.5

#Meeus formula 21.3
terrestrial_obliquity_coefficients = (
	s2d(23,26,21.448),
	-s2d(0,0,4680.93),
	-s2d(0,0,1.55),
	s2d(0,0,1999.25),
	-s2d(0,0,51.38),
	-s2d(0,0,249.67),
	-s2d(0,0,39.05),
	s2d(0,0,7.12),
	s2d(0,0,27.87),
	s2d(0,0,5.79),
	s2d(0,0,2.45)
)

#Adjust these coefficients for parameter T rather than U
terrestrial_obliquity_coefficients = [
	c * (1e-2) ** i for i,c in enumerate(terrestrial_obliquity_coefficients)
]

#Not entirely sure about this interpretation, but this is the difference
#between Meeus formulae 24.2 and 24.3 and seems to work
solar_perigee_coefficients = (
	280.46645 - 357.52910,
	36000.76932 - 35999.05030,
	0.0003032 + 0.0001559,
	0.00000048
)

#Meeus formula 24.2
solar_longitude_coefficients = (
	280.46645,
	36000.76983,
	0.0003032
)

#This value is taken from JPL Horizon and is essentially constant
lunar_inclination_coefficients = (
	5.145,
)

#Meeus formula 45.1
lunar_longitude_coefficients = (
	218.3164591,
	481267.88134236,
	-0.0013268,
	1/538841.0
	-1/65194000.0
)

#Meeus formula 45.7
lunar_node_coefficients = (
	125.0445550,
	-1934.1361849,
	0.0020762,
	1/467410.0,
	-1/60616000.0
)

#Meeus, unnumbered formula directly preceded by 45.7
lunar_perigee_coefficients = (
	83.3532430,
	4069.0137111,
	-0.0103238,
	-1/80053.0,
	1/18999000.0
)

#Now follow some useful auxiliary values, we won't need their speed.
#See notes on Table 6 in Schureman for I, nu, xi, nu', 2nu''
def _I(N, i, omega):
	N, i, omega = d2r * N, d2r*i, d2r*omega
	cosI = np.cos(i)*np.cos(omega)-np.sin(i)*np.sin(omega)*np.cos(N)
	return r2d*np.arccos(cosI)

def _xi(N, i, omega):
	N, i, omega = d2r * N, d2r*i, d2r*omega
	e1 = np.cos(0.5*(omega-i))/np.cos(0.5*(omega+i)) * np.tan(0.5*N)
	e2 = np.sin(0.5*(omega-i))/np.sin(0.5*(omega+i)) * np.tan(0.5*N)
	e1, e2 = np.arctan(e1), np.arctan(e2)
	e1, e2 = e1 - 0.5*N, e2 - 0.5*N
	return -(e1 + e2)*r2d

def _nu(N, i, omega):
	N, i, omega = d2r * N, d2r*i, d2r*omega
	e1 = np.cos(0.5*(omega-i))/np.cos(0.5*(omega+i)) * np.tan(0.5*N)
	e2 = np.sin(0.5*(omega-i))/np.sin(0.5*(omega+i)) * np.tan(0.5*N)
	e1, e2 = np.arctan(e1), np.arctan(e2)
	e1, e2 = e1 - 0.5*N, e2 - 0.5*N
	return (e1 - e2)*r2d

#Schureman equation 224
#Can we be more precise than B "the solar coefficient" = 0.1681?
def _nup(N, i, omega):
	I = d2r * _I(N, i, omega)
	nu = d2r * _nu(N, i, omega)
	return r2d * np.arctan(np.sin(2*I)*np.sin(nu)/(np.sin(2*I)*np.cos(nu)+0.3347))

#Schureman equation 232
def _nupp(N, i, omega):
	I = d2r * _I(N, i, omega)
	nu = d2r * _nu(N, i, omega)
	tan2nupp = (np.sin(I)**2*np.sin(2*nu))/(np.sin(I)**2*np.cos(2*nu)+0.0727)
	return r2d * 0.5 * np.arctan(tan2nupp)

AstronomicalParameter = namedtuple('AstronomicalParameter', ['value', 'speed'])

def astro(t):
	a = {}
	#We can use polynomial fits from Meeus to obtain good approximations to
	#some astronomical values (and therefore speeds).
	polynomials = {
			's':     lunar_longitude_coefficients,
			'h':     solar_longitude_coefficients,
			'p':     lunar_perigee_coefficients,
			'N':     lunar_node_coefficients,
			'pp':    solar_perigee_coefficients,
			'90':    (90.0,),
			'omega': terrestrial_obliquity_coefficients,
			'i':     lunar_inclination_coefficients
	}
	#Polynomials are in T, that is Julian Centuries; we want our speeds to be
	#in the more convenient unit of degrees per hour.
	dT_dHour = 1 / (24 * 365.25 * 100)
	for name, coefficients in polynomials.items():
		a[name] = AstronomicalParameter(
				np.mod(polynomial(coefficients, T(t)), 360.0),
				d_polynomial(coefficients, T(t)) * dT_dHour
		)

	#Some other parameters defined by Schureman which are dependent on the
	#parameters N, i, omega for use in node factor calculations. We don't need
	#their speeds.
	args = list(each.value for each in [a['N'], a['i'], a['omega']])
	for name, function in {
		'I':    _I,
		'xi':   _xi,
		'nu':   _nu,
		'nup':  _nup,
		'nupp': _nupp
	}.items():
		a[name] = AstronomicalParameter(np.mod(function(*args), 360.0), None)

	#We don't work directly with the T (hours) parameter, instead our spanning
	#set for equilibrium arguments #is given by T+h-s, s, h, p, N, pp, 90.
	#This is in line with convention.
	hour = AstronomicalParameter((JD(t) - np.floor(JD(t))) * 360.0, 15.0)
	a['T+h-s'] = AstronomicalParameter(
		hour.value + a['h'].value - a['s'].value,
		hour.speed + a['h'].speed - a['s'].speed
	)
	#It is convenient to calculate Schureman's P here since several node
	#factors need it, although it could be argued that these
	#(along with I, xi, nu etc) belong somewhere else.
	a['P'] = AstronomicalParameter(
		np.mod(a['p'].value -a['xi'].value,360.0),
		None
	)
	return a

