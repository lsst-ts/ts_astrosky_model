GH_PAGES_SOURCES = Makefile doc python README.rst HISTORY.rst tests SConstruct ups
SCONS_STUFF = config.log .sconsign.dblite .sconf_temp
BRANCH := $(shell git branch | grep \* | cut -d ' ' -f2)

.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

lint:
	flake8 python/lsst tests

test:
	py.test -v

test-all:
	tox

coverage:
	coverage run --source python/lsst -m unittest discover tests
	coverage report -m
	coverage html
	#open htmlcov/index.html

docs:
	rm -f doc/api/lsst*.rst
	rm -f doc/api/modules.rst
	sphinx-apidoc -H ts_astrosky_model -e -o doc/api python/lsst
	$(MAKE) -C doc clean
	$(MAKE) -C doc html
	#open doc/_build/html/index.html
