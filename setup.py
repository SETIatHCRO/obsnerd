#! /usr/bin/env python
# -*- mode: python; coding: utf-8 -*-
# Copyright 2024 David DeBoer
# Licensed under the 2-clause BSD license.

from setuptools import setup
import glob

setup_args = {
    'name': "obsnerd",
    'description': "NRDZ Observing for the ATA",
    'license': "BSD",
    'author': "David DeBoer",
    'author_email': "david.r.deboer@gmail.edu",
    'version': '0.3.0',
    'scripts': glob.glob('scripts/*'),
    'packages': ['obsnerd'],
    'include_package_data': True,
    'package_data': {"obsnerd": ['data/*.yaml']}
}

if __name__ == '__main__':
    setup(**setup_args)
