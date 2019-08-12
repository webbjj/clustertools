#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = open('README').read()
history = open('HISTORY').read().replace('.. :changelog:', '')

requirements = ['numpy', 'scipy>=0.14.0', 'galpy', 'numba', 'limepy']

setup(
    name='nbodypy',
    version='0.1.0',
    description='Package to Analyse N-body Simulations of Star Clusters',
    long_description=readme + '\n\n' + history,
    author='Jeremy Webb',
    author_email='webb@astro.utoronto.ca',
    url='https://github.com/webbjj/nbodypy.git',
    packages=[
        'nbodypy'
    ],
    package_dir={'main':'main','util':'util','custom':'custom'},
    license="BSD",
    keywords='nbodypy'
)
