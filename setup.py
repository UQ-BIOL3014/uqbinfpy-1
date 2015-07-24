#!/usr/bin/env python

import os
import sys
import glob

# Try and import pip. We'll stop if it is not present
try:
    import pip
except ImportError:
    print "Installation of uqbinfpy requires pip. Please install it!"
    print "See: http://pip.readthedocs.org/en/latest/installing.html"
    sys.exit(1)

from setuptools import setup

__title__ = 'uqbinfpy'
__version__ = '1.0.1'
__description__ = "uqbinfpy is a set of python modules for UQ BIOL3014"
__author__ = 'Mikael Boden'
__author_email__ = "m.boden@uq.edu.au"
__url__ = 'https://github.com/UQ-BIOL3014/uqbinfpy'


# Helper functions
if sys.argv[-1] == 'publish':
    print "Please use twine or do_release.sh"
    sys.exit()

if sys.argv[-1] == 'clean':
    os.system('rm -rf uqbinfpy.egg-info build dist')
    sys.exit()

packages = [__title__, ]

requires = []
with open('requirements.txt') as fin:
    lines = fin.readlines()
    for line in lines:
        requires.append(line.strip())

setup(
    name=__title__,
    version=__version__,
    description=__description__,
    long_description=open('README.rst').read(),
    author=__author__,
    author_email=__author_email__,
    url=__url__,
    packages=packages,
    test_suite="tests",
    package_dir={__title__: __title__},
    package_data={},
    #scripts=[__title__+'/TODO'],
    data_files=[('', ['requirements.txt', 'README.rst']),],
    include_package_data=True,
    install_requires=requires,
    zip_safe=False,
    classifiers=('Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Programming Language :: Python',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',),
)
