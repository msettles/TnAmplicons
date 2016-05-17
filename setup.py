#!/usr/bin/env python

'''
TnAmplicons
Analysis of Tn-Seq data, Transposon insertion site detection, initial
version is to process the samples (trim) primers (transoposon sequence)
detect the TA genomic insertion site and map the resulting files to the
genome.

Later version will exand on analysis
'''

import sys

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension


editdist = Extension('editdist', sources=['lib/editdist.c'])
trim = Extension('trim', sources=['lib/trim.c'])

try:
    version_num = open("VERSION", "r+").readline().strip()
except:
    sys.stderr.write("Error retrieving version_number")


config = \
    {
        'description': 'Processing of Illumina amplicon projects - TnSeq version',
        'author': 'Matt Settles',
        'url': 'https://github.com/msettles/TnAmplicons',
        'download_url': 'https://github.com/msettles/TnAmplicons',
        'author_email': 'settles@ucdavis.edu',
        'version': version_num,
        'install_requires': [],
        'packages': ['TnAmplicons'],
        'scripts': ['bin/TnAmplicons'],
        'name': 'TnAmplicons',
        "ext_package": 'TnAmplicons',
        'ext_modules': [editdist, trim]
    }

setup(**config)
