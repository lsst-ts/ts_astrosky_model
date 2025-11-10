# This file is part of ts_astrosky_model.
#
# Developed for the Vera Rubin Observatory Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

import math
import unittest
import warnings

import numpy
from astropy.time import Time
from lsst.ts.astrosky.model import AstronomicalSkyModel
from lsst.ts.dateloc import ObservatoryLocation
from rubin_scheduler.utils import survey_start_mjd


class AstronomicalSkyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Pick a time we know there should be skybrightness_pre data
        start_time = Time(survey_start_mjd(), format="mjd")
        cls.sky_start_stamp = start_time.unix
        cls.fixed_start_stamp = 1704582400
        warnings.filterwarnings("ignore", category=RuntimeWarning, append=True)
        warnings.filterwarnings("ignore", category=FutureWarning, append=True)
        cls.lsst_site = ObservatoryLocation()
        cls.lsst_site.for_lsst()
        cls.astro_sky = AstronomicalSkyModel(cls.lsst_site)
        cls.time_tolerance = 1e-6
        cls.sun_altitude = -12.0

    def create_ra_dec(self):
        self.ra_rads = numpy.radians(numpy.linspace(0.0, 90.0, 19))
        self.dec_rads = numpy.radians(numpy.linspace(-90.0, 0.0, 19))
        self.field_ids = numpy.arange(1, 20)

    def check_night_boundary_tuple(self, truth_set_timestamp, truth_rise_timestamp):
        (set_timestamp, rise_timestamp) = self.astro_sky.get_night_boundaries(
            self.sun_altitude
        )
        self.assertAlmostEqual(
            set_timestamp, truth_set_timestamp, delta=self.time_tolerance
        )
        self.assertAlmostEqual(
            rise_timestamp, truth_rise_timestamp, delta=self.time_tolerance
        )

    def test_basic_information_after_initial_creation(self):
        self.assertIsNotNone(self.astro_sky.date_profile)
        self.assertIsNotNone(self.astro_sky.sky_brightness_pre)
        self.assertIsNotNone(self.astro_sky.sun)
        self.assertTrue(self.astro_sky.exclude_planets)

    def test_update_mechanism(self):
        self.astro_sky.update(self.sky_start_stamp)
        self.assertEqual(self.astro_sky.date_profile.timestamp, self.sky_start_stamp)

    def test_sky_brightness_retrieval_internal_time_array_of_positions(self):
        self.create_ra_dec()
        self.astro_sky.update(self.sky_start_stamp)
        sky_mags = self.astro_sky.get_sky_brightness(self.ra_rads, self.dec_rads)
        self.assertEqual(len(sky_mags), 6)
        self.assertEqual(sky_mags["g"].size, self.field_ids.size)

    def test_sky_brightness_retrieval_from_timestamp_set_and_array_of_positions(self):
        initial_timestamp = self.sky_start_stamp
        time_step = 5.0 * 60.0
        number_of_steps = 10
        self.create_ra_dec()
        sky_mags = self.astro_sky.get_sky_brightness_timeblock(
            initial_timestamp, time_step, number_of_steps, self.ra_rads, self.dec_rads
        )
        self.assertEqual(len(sky_mags), number_of_steps)
        self.assertEqual(sky_mags[0]["g"].size, self.field_ids.size)

    def test_get_night_boundaries_20220101_sunset(self):
        # These night boundaries have hard-coded values based on an older
        # timestamp.

        # 2022/01/01
        # At sunset
        self.astro_sky.update(1641084532.843324)
        self.check_night_boundary_tuple(1641084532.89148, 1641113113.712405755558)

    def test_get_night_boundaries_20220101_night(self):
        # In night
        self.astro_sky.update(1641098823.29944)
        self.check_night_boundary_tuple(1641084532.89148, 1641113113.712405)

    def test_get_night_boundaries_20220102_sunrise(self):
        # At sunrise, next night bounds
        self.astro_sky.update(1641113113.755558)
        self.check_night_boundary_tuple(1641170940.944497, 1641199562.908091)

    def test_get_night_boundaries_20220102_daytime(self):
        # In daytime, next night bounds
        self.astro_sky.update(1641113114.755558)
        self.check_night_boundary_tuple(1641170940.944497, 1641199562.908091)
        self.astro_sky.update(1641133114.755558)
        self.check_night_boundary_tuple(1641170940.944497, 1641199562.908091)

    def test_get_night_boundaries_20220201(self):
        # 2022/02/01
        self.astro_sky.update(1643762299.348505)
        self.check_night_boundary_tuple(1643762299.384562, 1643793352.526473)

    def test_get_night_boundaries_20220208(self):
        # 2022/03/08
        self.astro_sky.update(1646784061.294245)
        self.check_night_boundary_tuple(1646784061.310997, 1646819228.773243)

    def test_get_night_boundaries_20220702(self):
        # 2022/07/02
        # At sunset
        self.astro_sky.update(1656802219.515093)
        self.check_night_boundary_tuple(1656802219.494713, 1656845034.721956)

    def test_get_night_boundaries_20220703(self):
        # At sunrise, next night bounds
        self.astro_sky.update(1656845034.721956)
        self.check_night_boundary_tuple(1656888641.705698, 1656931433.413172)

    def test_get_night_boundaries_20220703_daytime(self):
        # In daytime, next night bounds
        self.astro_sky.update(1656845035.696892)
        self.check_night_boundary_tuple(1656888641.705698, 1656931433.413172)

    def test_get_night_boundaries_20221017(self):
        # 2022/10/17
        self.astro_sky.update(1666050479.261601)
        self.check_night_boundary_tuple(1666050479.28522, 1666084046.849991)

    def test_get_night_boundaries_20250401(self):
        # 2025/04/01
        self.astro_sky.update(1743550264.401366)
        self.check_night_boundary_tuple(1743550264.405308, 1743588178.16701)

    def test_get_night_boundaries_20270621(self):
        # 2027/06/21
        self.astro_sky.update(1813618020.702736)
        self.check_night_boundary_tuple(1813618020.681721, 1813660970.015261)

    def test_get_night_boundaries_20310920(self):
        # 2031/09/20
        self.astro_sky.update(1947713387.331446)
        self.check_night_boundary_tuple(1947713387.340452, 1947750106.800034)

    def test_separation_function(self):
        initial_timestamp = self.fixed_start_stamp + (0.04166666666 * 3600 * 24)
        self.create_ra_dec()
        self.astro_sky.update(initial_timestamp)
        field_moon_sep = self.astro_sky.get_separation(
            "moon", self.ra_rads, self.dec_rads
        )
        self.assertEqual(field_moon_sep.size, 19)
        # self.assertAlmostEqual(field_moon_sep[0], numpy.radians(64.6988849))
        field_sun_sep = self.astro_sky.get_separation(
            "sun", self.ra_rads, self.dec_rads
        )
        self.assertEqual(field_sun_sep.size, 19)
        # self.assertAlmostEqual(field_sun_sep[0], numpy.radians(67.06949045))

    def test_moon_sun_information(self):
        initial_timestamp = self.fixed_start_stamp + (0.04166666666 * 3600 * 24)
        self.create_ra_dec()
        self.astro_sky.update(initial_timestamp)
        info = self.astro_sky.get_moon_sun_info(self.ra_rads, self.dec_rads)
        self.assertEqual(len(info), 11)
        self.assertAlmostEqual(info["moonPhase"], 32.01068613858965, delta=1e-6)
        self.assertEqual(len(info["moonDist"]), self.ra_rads.size)
        self.assertAlmostEqual(info["moonDist"][0], 1.2303584595787405, delta=1e-7)
        self.assertAlmostEqual(info["moonDec"], -0.3404378672161562, delta=1e-7)
        self.assertAlmostEqual(info["moonRA"], 3.9285205954901437, delta=1e-7)
        self.assertEqual(len(info["solarElong"]), self.ra_rads.size)
        self.assertAlmostEqual(info["solarElong"][0], 1.1781286433066063, delta=1e-7)

    def test_target_information(self):
        initial_timestamp = self.fixed_start_stamp + (0.04166666666 * 3600 * 24)
        self.create_ra_dec()
        self.astro_sky.update(initial_timestamp)
        info = self.astro_sky.get_target_information(self.ra_rads, self.dec_rads)
        self.assertEqual(len(info), 3)
        self.assertEqual(info["airmass"].size, self.field_ids.size)
        self.assertEqual(info["altitude"].size, self.ra_rads.size)
        self.assertEqual(info["azimuth"].size, self.ra_rads.size)
        self.assertAlmostEqual(info["airmass"][0], 1.9853363754043085, delta=1e-5)
        self.assertAlmostEqual(info["altitude"][0], 0.52786436029017303, delta=1e-5)
        self.assertFalse(numpy.isnan(info["azimuth"][0]))

    def test_configure(self):
        self.astro_sky.configure(exclude_planets=False)
        self.assertFalse(self.astro_sky.exclude_planets)

    def test_update_location(self):
        # Kitt Peak
        latitude = 31.959954
        longitude = -111.599808
        height = 2097.0
        location = ObservatoryLocation(
            math.radians(latitude), math.radians(longitude), height
        )
        self.astro_sky.update_location(location)
        dp = self.astro_sky.date_profile
        self.assertEqual(dp.location.latitude, latitude)
        self.assertEqual(dp.location.longitude, longitude)
        self.assertEqual(dp.location.height, height)

    def test_get_hour_angle(self):
        initial_timestamp = self.fixed_start_stamp + (0.04166666666 * 3600 * 24)
        self.create_ra_dec()
        self.astro_sky.update(initial_timestamp)
        hour_angle = self.astro_sky.get_hour_angle(self.ra_rads)
        self.assertEqual(hour_angle.size, self.ra_rads.size)
        self.assertAlmostEqual(hour_angle[0], 0.6455681350636748, delta=1e-6)
