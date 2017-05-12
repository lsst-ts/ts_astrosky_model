============
Installation
============

The installation of ``ts_astrosky_model`` requires the use of the ``ts_dateloc`` package. Installation instructions for that package can be found `here <https://github.com/lsst-ts/ts_dateloc/blob/master/doc/installation.rst>`_. Once those instructions are complete, install the source code into your favorite location (called ``gitdir``) via::

	git clone https://github.com/lsst-ts/ts_astrosky_model.git

With the stack environment setup as instructed, declare the package to EUPS::

	cd gitdir/ts_astrosky_model
	eups declare ts_astrosky_model git -r . -c
	setup ts_astrosky_model git
	scons