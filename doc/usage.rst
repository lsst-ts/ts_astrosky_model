=====
Usage
=====

The main class for use in this package is the :py:class:`.AstronomicalSkyModel`. It is a wrapper for the ``sims_skybrightness_pre`` package, but provides other convenience methods. It requires a :py:class:`lsst.ts.dateloc.ObservatoryLocation` instance. The example will start with creating a LSST instance.

.. code-block:: python

  from lsst.ts.dateloc import ObservatoryLocation
  lsst = ObservatoryLocation()
  lsst.for_lsst()

Next, create an instance of the astronomical sky model.

.. code-block:: python

  from lsst.ts.astrosky.model import AstronomicalSkyModel
  asky = AstronomicalSkyModel(lsst)

The model needs to be updated to a given timestamp for calculations. This is handled by doing this.

.. code-block:: python

  asky.update(1672542000)

The astronomical sky model uses OpSim field Ids for gathering sky brightness information. The Ids range from 1 to 5292. Field Id = 1 contains the South Celestial Pole, so from the LSST site, it is always visible. Retrieve the sky brightness for that field. 

.. code-block:: python

  import numpy
  ids = numpy.array([1])
  sb = asky.get_sky_brightness(ids)
  sb
  {'g': array([ 19.5263785]),
   'i': array([ 19.24296218]),
   'r': array([ 19.47577689]),
   'u': array([ 20.17420233]),
   'y': array([ 17.52278367]),
   'z': array([ 18.59328021])}

The model can be used to find much more information. Using the SCP field, find the airmass of the field.

.. code-block:: python 

  airmass = asky.get_airmass(ids)
  airmass
  array([ 1.98534993])

Now, find the altitude, azimuth and hour angle of the SCP field.

.. code-block:: python

  ra = numpy.radians(numpy.array([0.0]))
  dec = numpy.radians(numpy.array([-90.0]))
  ha = asky.get_hour_angle(ra)
  array([ 1.30489771])
  alt, az = asky.get_alt_az(ra, dec)
  alt
  array([ 0.52786436])
  az
  array([ 3.14159265])

Now, find Moon and Sun information. This includes separation values from the field to the Moon and Sun.

.. code-block:: python

  info = asky.get_moon_sun_info(ra, dec)
  info
  {'moonAlt': array([ 0.53168999]),
   'moonAz': array([ 5.4229169]),
   'moonDec': 0.23239105825985337,
   'moonDist': array([ 1.80318739]),
   'moonPhase': 71.086377626099321,
   'moonRA': 0.56876998782850385,
   'solarElong': array([ 1.168883]),
   'sunAlt': array([-0.53997986]),
   'sunAz': array([ 3.64194876]),
   'sunDec': -0.4019133262390896,
   'sunRA': 4.9100206262887278}

See the API documentation for :py:class:`.AstronomicalSkyModel`.