#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2021, MeteoSwiss, the authors

from setuptools import setup, find_packages

requirements = [
    "netcdf4",
    "numba",
    "dask",
    "numpy",
    "pandas",
    "scipy",
    "utm",
    "xarray",
    "yaconfigobject",
]

setup_requirements = [
    "setuptools_scm",
]

test_requirements = [
    "pytest",
]

extras = {
    "test": test_requirements,
}

packages = find_packages(include=["topo_descriptors"])

package_dir = {}

package_data = {}

setup(
    name="topo-descriptors",
    packages=packages,
    use_scm_version=True,
    author="Mathieu Schaer",
    author_email="mathieu.schaer@meteoswiss.ch",
    maintainer="Daniele Nerini",
    maintainer_email="daniele.nerini@meteoswiss.ch",
    description="A library to compute DEM-based topographical descriptors.",
    long_description=open("README.md").read() + "\n\n" + open("HISTORY.rst").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/MeteoSwiss/topo-descriptors",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
        "Operating System :: OS Independent",
    ],
    license="BSD-3-Clause license",
    keywords="topo_descriptors",
    entry_points={},
    py_modules=["topo-descriptors"],
    include_package_data=True,
    install_requires=requirements,
    package_dir=package_dir,
    package_data=package_data,
    setup_requires=setup_requirements,
    tests_require=test_requirements,
    extras_require=extras,
)
