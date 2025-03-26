.. :changelog:

History
-------

1.4.1 (2025-03-26)
~~~~~~~~~~~~~~~~~~

* Update to work with astropy 7.


1.4.0 (2023-01-26)
~~~~~~~~~~~~~~~~~~

* Upgrade to work with rubin-sim >1.

1.3.0 (2022-06-23)
~~~~~~~~~~~~~~~~~~

* Upgrade to pyproject.toml.
* Update unit tests to be compatible with the latest version of rubin-sim.
  Basically, had to change the time when the tests are exe
* Update AstroSky to make it compatible with the latest version of rubin-sim.

1.2.0 (2021-06-22)
~~~~~~~~~~~~~~~~~~

* Remove unused `__init__` files.
* Add license header to all files.
* Reformat code with black and add latest style configurations.
* Replace lsst.sims with rubin_sim and remove import of future division.
* Add conda packaging.

1.1.1 (2018-11-02)
~~~~~~~~~~~~~~~~~~

* Update unit tests.

1.1.0 (2018-04-19)
~~~~~~~~~~~~~~~~~~

* Implement ra/dec on get_sky_brightness_timeblock and fix unit tests.
* Cleaning up dependency on OpSim field ID.

1.0.0 (2017-05-22)
~~~~~~~~~~~~~~~~~~

* Provides access to astronomical sky related information:

  * sky brightness
  * airmass
  * sun and moon information
  * target information

    * altitude/azimuth
    * moon separation
    * sun separation

  * hour angle
  * night boundaries
