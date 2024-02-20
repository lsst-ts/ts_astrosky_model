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

__all__ = ["AstronomicalSkyModel"]

import logging
import typing
import warnings

import numpy
import palpy
from astropy import units
from astropy.coordinates import get_moon, get_sun
from astropy.time import Time
from lsst.ts.astrosky.model.sun import Sun
from lsst.ts.dateloc import DateProfile
from rubin_scheduler import __version__ as sky_model_pre_version
from rubin_scheduler.skybrightness_pre import SkyModelPre
from rubin_scheduler.utils import _ra_dec2_hpid


class AstronomicalSkyModel(object):
    def __init__(self, location):
        """Initialize the class.

        Parameters
        ----------
        location : `lsst.ts.dateloc.ObservatoryLocation`
            The instance containing the observatory location information.
        """
        self.log = logging.getLogger("sky_model.AstronomicalSkyModel")
        self.date_profile = DateProfile(0, location)
        self.sky_brightness_pre = SkyModelPre()
        self._sb_nside = self.sky_brightness_pre.nside
        self.sun = Sun()
        self.exclude_planets = True

    def configure(self, exclude_planets):
        """Add configuration for the sky brightness model.

        Parameters
        ----------
        exclude_planets : `bool`
            Flag to mask planets in sky brightness information.
        """
        self.exclude_planets = exclude_planets

    def get_airmass(self, ra, dec):
        """Get the airmass of the fields.

        Parameters
        ----------
        ra : `numpy.array`
            The right ascension (radians) of the sky position.
        dec : `numpy.array`
            The declination (radians) of the sky position.

        Returns
        -------
        `numpy.array`
            The set of airmasses.
        """
        alts, az = self.get_alt_az(ra, dec)
        # Calculate airmass as in the ESO model / rubin_sim.Skybrightness
        airmass = 1.0 / numpy.cos(numpy.pi / 2.0 - alts)
        return airmass

    def get_alt_az(self, ra, dec):
        """Get the altitude (radians) and azimuth (radians) of a given sky
        position.

        Parameters
        ----------
        ra : `numpy.array`
            The right ascension (radians) of the sky position.
        dec : `numpy.array`
            The declination (radians) of the sky position.

        Returns
        -------
        `tuple`
            The altitude and azimuth of the sky position.
        """
        hour_angle = self.date_profile.lst_rad - ra
        azimuth, altitude = palpy.de2hVector(
            hour_angle, dec, self.date_profile.location.latitude_rad
        )
        return altitude, azimuth

    def get_hour_angle(self, ra):
        """Get the hour angle (radians) of a sky position.

        Parameters
        ----------
        ra : `numpy.array`
            The right ascension (radians) of the sky position.

        Returns
        -------
        `numpy.array`
            The hour angle (radians) of the sky position on range of -pi to pi.
        """
        hour_angle = self.date_profile.lst_rad - ra
        hour_angle = numpy.where(
            hour_angle < -numpy.pi, hour_angle + (2.0 * numpy.pi), hour_angle
        )
        hour_angle = numpy.where(
            hour_angle > numpy.pi, hour_angle - (2.0 * numpy.pi), hour_angle
        )
        return hour_angle

    def get_moon_sun_info(self, field_ra, field_dec):
        """Return the current moon and sun information.

        This function gets the right ascension, declination, altitude, azimuth,
        phase and angular distance from target (by given ra and dec) for the
        moon and the right ascension, declination, altitude, azimuth and solar
        elongation from target (by given ra and dec) for the sun.

        Parameters
        ----------
        field_ra : `numpy.array`
            The target right ascension (radians).
        field_dec : `numpy.array`
            The target declination (radians).

        Returns
        -------
        `dict`
            The set of information pertaining to the moon and sun. All angles
            are in radians.
        """
        attrs = self.get_sun_moon()
        moon_distance = self.get_separation("moon", field_ra, field_dec)
        sun_distance = self.get_separation("sun", field_ra, field_dec)

        keys = ["moonRA", "moonDec", "sunRA", "sunDec"]
        info_dict = {}
        for key in keys:
            info_dict[key] = attrs[key]
        # moonSunSep is in degrees! Oops!
        info_dict["moonPhase"] = attrs["moonSunSep"] / 180.0 * 100.0
        info_dict["moonAlt"], info_dict["moonAz"] = self.get_alt_az(
            numpy.array([attrs["moonRA"]]), numpy.array([attrs["moonDec"]])
        )
        info_dict["sunAlt"], info_dict["sunAz"] = self.get_alt_az(
            numpy.array([attrs["sunRA"]]), numpy.array([attrs["sunDec"]])
        )
        info_dict["moonDist"] = moon_distance
        info_dict["solarElong"] = sun_distance
        return info_dict

    def get_night_boundaries(
        self, sun_altitude, upper_limb_correction=False, precision=6
    ):
        """Return the set/rise times of the sun for the given altitude.

        This function calculates the night boundaries (the set and rise times)
        for a given sun altitude. It uses the currently stored timestamp in the
        lsst.ts.dateloc.DateProfile instance.

        Parameters
        ----------
        sun_altitude : `float`
            The altitude of the sun to get the set/rise times for.
        upper_limb_correction : `bool`
            Set to True is the upper limb correction should be calculated.
        precision : `int`, optional
            The place to round the rise/set times.

        Returns
        -------
        `tuple` (`float`, `float`)
            A tuple of the set and rise times, respectively, for the sun
            altitude.
        """
        longitude, latitude = (
            self.date_profile.location.longitude,
            self.date_profile.location.latitude,
        )

        current_timestamp = self.date_profile.timestamp

        midnight_timestamp = self.date_profile.midnight_timestamp()
        (rise_time, set_time) = self.sun.altitude_times(
            midnight_timestamp, longitude, latitude, sun_altitude, upper_limb_correction
        )

        set_timestamp = midnight_timestamp + (
            set_time * self.date_profile.SECONDS_IN_HOUR
        )
        rise_timestamp = midnight_timestamp + (
            rise_time * self.date_profile.SECONDS_IN_HOUR
        )

        if current_timestamp < rise_timestamp:
            midnight_timestamp = self.date_profile.previous_midnight_timestamp()
            (_, set_time) = self.sun.altitude_times(
                midnight_timestamp,
                longitude,
                latitude,
                sun_altitude,
                upper_limb_correction,
            )

            set_timestamp = midnight_timestamp + (
                set_time * self.date_profile.SECONDS_IN_HOUR
            )

        else:
            midnight_timestamp = self.date_profile.next_midnight_timestamp()
            (rise_time, _) = self.sun.altitude_times(
                midnight_timestamp,
                longitude,
                latitude,
                sun_altitude,
                upper_limb_correction,
            )

            rise_timestamp = midnight_timestamp + (
                rise_time * self.date_profile.SECONDS_IN_HOUR
            )

        return (round(set_timestamp, precision), round(rise_timestamp, precision))

    def get_separation(self, body, field_ra, field_dec):
        """Return the separation between a body and a set of field coordinates.

        This function returns the separation (in radians) between the given
        body (either moon or sun) and a given set of fields. It uses a list of
        (RA, Dec) coordinates. This function assumes that
        :meth:`.get_sky_brightness` has been run.

        Parameters
        ----------
        body : `str`
            The name of the body to calculate the separation. Either moon or
            sun.
        field_ra : `numpy.array` of `float`
            The list of field right ascensions in radians.
        field_dec : `numpy.array` of `float`
            The list of field declinations in radians.

        Returns
        -------
        `numpy.array` of `float`
            The list of field-moon separations in radians.
        """
        attrs = self.get_sun_moon()
        return palpy.dsepVector(
            field_ra,
            field_dec,
            numpy.full_like(field_ra, attrs["{}RA".format(body)]),
            numpy.full_like(field_dec, attrs["{}Dec".format(body)]),
        )

    def get_sky_brightness(
        self, ra, dec, extrapolate=False, override_exclude_planets=None
    ):
        """Get the LSST 6 filter sky brightness for a set of coordinates at a
        single time.

        This function retrieves the LSST 6 filter sky brightness magnitudes for
        a given set of ra/decs at the MJD kept by the
        lsst.ts.dateloc.DateProfile.

        Parameters
        ----------
        ra : `numpy.array`
            The set of fields RA coordinates to retrieve the sky brightness
            for.
        dec : `numpy.array`
            The set of fields RA coordinates to retrieve the sky brightness
            for.
        extrapolate : `boolean`, optional
            Flag to extrapolate fields with bad sky brightness to nearest field
            that is good.
        override_exclude_planets : `boolean`, optional
            (Deprecated) Override the internally stored exclude_planets flag.

        Returns
        -------
        `numpy.ndarray`
            The LSST 6 filter sky brightness magnitudes.
        """
        if override_exclude_planets is not None:
            warnings.warn(
                "Parameter `override_exclude_planets` is deprecated.",
                DeprecationWarning,
            )

        ids = _ra_dec2_hpid(self._sb_nside, ra, dec)

        return self.sky_brightness_pre.return_mags(
            self.date_profile.mjd,
            indx=ids,
            badval=float("nan"),
            extrapolate=extrapolate,
        )

    def get_sky_brightness_timeblock(self, timestamp, timestep, num_steps, ra, dec):
        """Get LSST 6 filter sky brightness for a set of fields for a range of
        times.

        This function retrieves the LSST 6 filter sky brightness magnitudes for
        a given set of ra/decs at a range of MJDs provided via the timeblock
        information.


        Parameters
        ----------
        timestamp : `float`
            The UTC timestamp to start the time block.
        timestep : `float`
            The number of seconds to increment the timestamp with.
        num_steps : `int`
            The number of steps to create for the time block.
        ra : `numpy.array`
            The set of fields RA coordinates to retrieve the sky brightness
            for.
        dec : `numpy.array`
            The set of fields RA coordinates to retrieve the sky brightness
            for.

        Returns
        -------
        `numpy.ndarray`
            The LSST 6 filter sky brightness magnitudes.
        """
        dp = DateProfile(0, self.date_profile.location)
        mags = []
        ids = _ra_dec2_hpid(self._sb_nside, ra, dec)

        for i in range(num_steps):
            ts = timestamp + i * timestep
            mjd, _ = dp(ts)
            mags.append(
                self.sky_brightness_pre.return_mags(
                    dp.mjd,
                    indx=ids,
                    badval=float("nan"),
                    extrapolate=False,
                )
            )

        return mags

    def get_target_information(self, ra, dec):
        """Get information about target(s).

        This function gathers altitude (radians) and azimuth (radians)
        information for the target.

        Parameters
        ----------
        ra : `numpy.array`
            The field right ascension (radians).
        dec : `numpy.array`
            The field declination (radians).

        Returns
        -------
        `dict`
            Set of information about the target(s).
        """
        info_dict = {}
        altitude, azimuth = self.get_alt_az(ra, dec)
        info_dict["airmass"] = 1.0 / numpy.cos(numpy.pi / 2.0 - altitude)
        info_dict["altitude"] = altitude
        info_dict["azimuth"] = azimuth
        return info_dict

    def sky_brightness_config(self):
        """Get the configuration from the SkyModelPre files.

        Returns
        -------
        `list` of `tuple` (key, value)
        """
        config = []
        header = self.sky_brightness_pre.header
        config.append(("sky_brightness_pre/program_version", sky_model_pre_version))
        config.append(("sky_brightness_pre/file_version", header["version"]))
        config.append(("sky_brightness_pre/fingerprint", header["fingerprint"]))
        config.append(("sky_brightness_pre/moon_dist_limit", header["moon_dist_limit"]))
        config.append(
            ("sky_brightness_pre/planet_dist_limit", header["planet_dist_limit"])
        )
        config.append(("sky_brightness_pre/airmass_limit", header["airmass_limit"]))
        config.append(("sky_brightness_pre/timestep", header["timestep"] * 24 * 3600))
        config.append(
            ("sky_brightness_pre/timestep_max", header["timestep_max"] * 24 * 3600)
        )
        config.append(("sky_brightness_pre/dm", header["dm"]))
        return config

    def update(self, timestamp):
        """Update the internal timestamp.

        Parameters
        ----------
        timestamp : `float`
            The UTC timestamp to update the internal timestamp to.
        """
        self.date_profile.update(timestamp)

    def update_location(self, location):
        """Update the current location.

        Parameters
        ----------
        location : `lsst.ts.dateloc.ObservatoryLocation`
            The instance of an ObservatoryLocation containing the relevant
            information.
        """
        self.date_profile.location = location

    def get_sun_moon(self) -> typing.Dict[str, float]:
        """Return sun and moon position.

        Returns
        -------
        `dict`
            Dictionary with the sun and moon position.
        """
        mjd = Time(self.date_profile.mjd, format="mjd")

        coord_sun = get_sun(mjd)
        coord_moon = get_moon(mjd)
        moon_sun_sep = numpy.degrees(
            palpy.dsepVector(
                numpy.array([coord_sun.ra.to(units.rad).value]),
                numpy.array([coord_sun.dec.to(units.rad).value]),
                numpy.array([coord_moon.ra.to(units.rad).value]),
                numpy.array([coord_moon.dec.to(units.rad).value]),
            )
        )[0]
        return dict(
            sunRA=coord_sun.ra.to(units.rad).value,
            sunDec=coord_sun.dec.to(units.rad).value,
            moonRA=coord_moon.ra.to(units.rad).value,
            moonDec=coord_moon.dec.to(units.rad).value,
            moonSunSep=moon_sun_sep,
        )
