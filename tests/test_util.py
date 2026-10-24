#!/usr/bin/env python

import datetime
import os.path
import shutil
import tempfile
import pytest
from urllib.error import URLError

import numpy as np

import clawpack.geoclaw.util as util
from clawpack.geoclaw.util import NOAA_API_URL, fetch_noaa_tide_data

@pytest.mark.python
class TestFetchNoaaTideData:
    station = '1234567'
    begin_date = datetime.datetime(2000, 10, 30, 12, 0)
    end_date = datetime.datetime(2000, 10, 30, 12, 24)

    water_level_url = '{}?product={}&station={}'.format(
        NOAA_API_URL, 'water_level', station)
    predictions_url = '{}?product={}&station={}'.format(
        NOAA_API_URL, 'predictions', station)

    @classmethod
    def setup_class(cls):
        cls.temp_dir = tempfile.mkdtemp()

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.temp_dir)

    def test_retrieve_and_cache(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_retrieve_and_cache.__name__)

        water_level_response = \
            ('Date Time, Water Level, Sigma, O or I (for verified), F, R, L, Quality\n'
             '2000-10-30 12:00,1.001,0.001,0,0,0,0,v\n'
             '2000-10-30 12:06,1.002,0.002,0,0,0,0,v\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        predictions_response = ('Date Time, Prediction\n'
                                '2000-10-30 12:00,1.101\n'
                                '2000-10-30 12:06,1.102\n'
                                '2000-10-30 12:12,1.103\n'
                                '2000-10-30 12:18,1.104\n'
                                '2000-10-30 12:24,1.105\n')

        date_time_range = ['2000-10-30T12:00',
                           '2000-10-30T12:06',
                           '2000-10-30T12:12',
                           '2000-10-30T12:18',
                           '2000-10-30T12:24']

        date_time_expected = np.array(date_time_range, dtype='datetime64[m]')
        water_level_expected = np.array([1.001, 1.002, 1.003, 1.004, 1.005])
        prediction_expected = np.array([1.101, 1.102, 1.103, 1.104, 1.105])

        # monkey patch urllib to return mock data
        def mock_read_response(url):
            if 'product=water_level' in url:
                return water_level_response
            if 'product=predictions' in url:
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        # first time, should fetch data and save in cache
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir)

        # configure endpoint to not return data anymore
        self._monkey_patch_urlopen(lambda url: 'Should read from cache')

        # second time, should read data from cache
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir)

    def test_api_error(self):
        cache_dir = os.path.join(self.temp_dir, self.test_api_error.__name__)

        # configure endpoint to return an error
        self._monkey_patch_urlopen(lambda url: 'Something went wrong')

        # should return None
        d,w,p = fetch_noaa_tide_data(self.station, self.begin_date, self.end_date,
                             cache_dir=cache_dir)
        assert d == None, '*** expected d == None'

    def test_date_time_range_mismatch(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_date_time_range_mismatch.__name__)

        # missing first two entries
        water_level_response = \
            ('Date Time, Water Level, Sigma, O or I (for verified), F, R, L, Quality\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        # missing last two entries
        predictions_response = ('Date Time, Prediction\n'
                                '2000-10-30 12:00,1.101\n'
                                '2000-10-30 12:06,1.102\n'
                                '2000-10-30 12:12,1.103\n')

        # monkey patch urllib to return mock data
        def mock_read_response(url):
            if 'product=water_level' in url:
                return water_level_response
            if 'product=predictions' in url:
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        # should raise ValueError
        with pytest.raises(ValueError):
            fetch_noaa_tide_data(self.station, self.begin_date, self.end_date,
                                 cache_dir=cache_dir)

    def test_missing_values(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_missing_values.__name__)

        # missing fist two water level values
        water_level_response = \
            ('Date Time, Water Level, Sigma, O or I (for verified), F, R, L, Quality\n'
             '2000-10-30 12:00,,0.001,0,0,0,0,v\n'
             '2000-10-30 12:06,,0.002,0,0,0,0,v\n'
             '2000-10-30 12:12,1.003,0.003,0,0,0,0,v\n'
             '2000-10-30 12:18,1.004,0.004,0,0,0,0,v\n'
             '2000-10-30 12:24,1.005,0.005,0,0,0,0,v\n')

        # missing a prediction value
        predictions_response = ('Date Time, Prediction\n'
                                '2000-10-30 12:00,1.101\n'
                                '2000-10-30 12:06,1.102\n'
                                '2000-10-30 12:12,1.103\n'
                                '2000-10-30 12:18,\n'
                                '2000-10-30 12:24,1.105\n')

        date_time_range = ['2000-10-30T12:00',
                           '2000-10-30T12:06',
                           '2000-10-30T12:12',
                           '2000-10-30T12:18',
                           '2000-10-30T12:24']

        date_time_expected = np.array(date_time_range, dtype='datetime64[m]')
        water_level_expected = np.array([np.nan, np.nan, 1.003, 1.004, 1.005])
        prediction_expected = np.array([1.101, 1.102, 1.103, np.nan, 1.105])

        # monkey patch urllib to return mock data
        def mock_read_response(url):
            if 'product=water_level' in url:
                return water_level_response
            if 'product=predictions' in url:
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir)

    def test_hourly_height(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_hourly_height.__name__)

        # hourly heights (pre-1996 records) with a missing value; NOAA's
        # header carries a trailing space
        hourly_height_response = \
            ('Date Time, Water Level, Sigma, I, L \n'
             '1992-12-11 00:00,2.835,0.0,0,0\n'
             '1992-12-11 01:00,,0.0,0,0\n'
             '1992-12-11 02:00,3.214,0.0,0,0\n')

        predictions_response = ('Date Time, Prediction\n'
                                '1992-12-11 00:00,2.208\n'
                                '1992-12-11 01:00,2.331\n'
                                '1992-12-11 02:00,2.250\n')

        date_time_expected = np.array(['1992-12-11T00:00',
                                       '1992-12-11T01:00',
                                       '1992-12-11T02:00'],
                                      dtype='datetime64[m]')
        water_level_expected = np.array([2.835, np.nan, 3.214])
        prediction_expected = np.array([2.208, 2.331, 2.250])

        # predictions must be requested hourly to match the hourly heights
        def mock_read_response(url):
            if 'product=hourly_height' in url:
                return hourly_height_response
            if 'product=predictions' in url:
                assert 'interval=h' in url
                return predictions_response
            raise AssertionError
        self._monkey_patch_urlopen(mock_read_response)

        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir,
                               product='hourly_height')

        # hourly predictions are cached apart from 6 minute ones
        assert os.path.isdir(os.path.join(cache_dir, 'predictions_h'))
        assert not os.path.exists(os.path.join(cache_dir, 'predictions'))

        # second time, should read data from cache
        self._monkey_patch_urlopen(lambda url: 'Should read from cache')
        self._fetch_and_assert(date_time_expected, water_level_expected,
                               prediction_expected, cache_dir,
                               product='hourly_height')

    def test_network_failure(self):
        cache_dir = os.path.join(self.temp_dir,
                                 self.test_network_failure.__name__)

        # unreachable network (e.g. an HPC compute node): report, don't raise
        def mock_urlopen(url, timeout=None):
            assert timeout == 5
            raise URLError('Network is unreachable')
        util.urlopen = mock_urlopen

        d, w, p = fetch_noaa_tide_data(self.station, self.begin_date,
                                       self.end_date, cache_dir=cache_dir,
                                       product='hourly_height', timeout=5)
        assert d is None and w is None and p is None

    def test_unknown_product(self):
        with pytest.raises(ValueError):
            fetch_noaa_tide_data(self.station, self.begin_date, self.end_date,
                                 cache_dir=self.temp_dir, product='monthly_mean')

    def _fetch_and_assert(self, date_time_expected, water_level_expected,
                          prediction_expected, cache_dir,
                          product='water_level'):
        # fetch data
        date_time, water_level, prediction = fetch_noaa_tide_data(
            self.station, self.begin_date, self.end_date, cache_dir=cache_dir,
            product=product)

        # make sure data was correctly retrieved
        np.testing.assert_equal(date_time, date_time_expected)
        np.testing.assert_equal(water_level, water_level_expected)
        np.testing.assert_equal(prediction, prediction_expected)

    @staticmethod
    def _monkey_patch_urlopen(mock_read_response):
        """Modify 'util.urlopen' to return a custom response."""
        def mock_urlopen(url, timeout=None):
            class MockHttpResponse:
                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc_value, traceback):
                    pass

                def read(self):
                    # encode mock response
                    text = mock_read_response(url)
                    return text.encode('utf-8')

            return MockHttpResponse()

        # replace 'util.urlopen' with mock implementation
        util.urlopen = mock_urlopen
