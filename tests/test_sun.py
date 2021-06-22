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

from datetime import datetime
import unittest

from lsst.ts.astrosky.model import Sun

"""Set timestamp as 2022-01-01 0h UTC"""
LSST_START_TIMESTAMP = 1640995200.0
"""LSST latitude (degrees)"""
LSST_LATITUDE = -30.2444
"""LSST longitude (degrees)"""
LSST_LONGITUDE = -70.7494

class SunTest(unittest.TestCase):

    def setUp(self):
        self.sun = Sun()
        self.timestamp = LSST_START_TIMESTAMP
        self.date = datetime.utcfromtimestamp(self.timestamp)
        self.days = 8037.0
        self.latitude = LSST_LATITUDE
        self.longitude = LSST_LONGITUDE
        # Set an arbitrary sun altitude for rise/set calculations
        self.arb_alt = -15.0
        # Tolerance for float comparisions
        self.tol = 1E-8

    def test_days_since_2000_jan_0(self):
        days = self.sun.days_since_2000_jan_0(self.timestamp)
        self.assertEqual(days, self.days)

    def test_angle_normalizations(self):
        angle1 = 540.0
        self.assertEqual(self.sun.normalize(angle1), 180.0)
        self.assertEqual(self.sun.normalize(angle1, from_minus_180=True), -180.0)

    def test_gmst0(self):
        self.assertAlmostEqual(self.sun.gmst0(self.days), 100.63516802400045, delta=self.tol)

    def test_position(self):
        solar_lon, _ = self.sun.position(self.days)
        self.assertAlmostEqual(solar_lon, 280.5436722457089, delta=self.tol)

    def test_ra_dec(self):
        (ra, dec, _) = self.sun.ra_dec(self.days)
        self.assertAlmostEqual(ra, -78.53240325724865, delta=self.tol)
        self.assertAlmostEqual(dec, -23.017734721208154, delta=self.tol)

    def test_altitude_times_arbitrary_altitude_no_limb_correction(self):
        (sunrise, sunset) = self.sun.altitude_times(self.timestamp, self.longitude, self.latitude,
                                                    self.arb_alt, False)
        self.assertAlmostEqual(sunrise, 8.449821291129762, delta=self.tol)
        self.assertAlmostEqual(sunset, 25.105364046890113, delta=self.tol)

    def test_altitude_times_arbitrary_altitude_with_limb_correction(self):
        (sunrise, sunset) = self.sun.altitude_times(self.timestamp, self.longitude, self.latitude,
                                                    self.arb_alt, True)
        self.assertAlmostEqual(sunrise, 8.423005170958238, delta=self.tol)
        self.assertAlmostEqual(sunset, 25.132180167061637, delta=self.tol)

    def test_rise_set(self):
        (sunrise, sunset) = self.sun.rise_set(self.timestamp, self.longitude, self.latitude)
        self.assertAlmostEqual(sunrise, 9.750138553498452, delta=self.tol)
        self.assertAlmostEqual(sunset, 23.805046784521423, delta=self.tol)

    def test_civil_twilight(self):
        (sunrise, sunset) = self.sun.civil_twilight(self.timestamp, self.longitude, self.latitude)
        self.assertAlmostEqual(sunrise, 9.29521962888819, delta=self.tol)
        self.assertAlmostEqual(sunset, 24.259965709131684, delta=self.tol)

    def test_nautical_twilight(self):
        (sunrise, sunset) = self.sun.nautical_twilight(self.timestamp, self.longitude, self.latitude)
        self.assertAlmostEqual(sunrise, 8.740506636842959, delta=self.tol)
        self.assertAlmostEqual(sunset, 24.814678701176916, delta=self.tol)

    def test_astronomical_twilight(self):
        (sunrise, sunset) = self.sun.astronomical_twilight(self.timestamp, self.longitude, self.latitude)
        self.assertAlmostEqual(sunrise, 8.14717205365863, delta=self.tol)
        self.assertAlmostEqual(sunset, 25.408013284361246, delta=self.tol)
