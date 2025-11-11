import setuptools_scm
from setuptools import setup

setup(
    version=setuptools_scm.get_version(
        write_to="python/lsst/ts/astrosky/model/version.py"
    )
)
