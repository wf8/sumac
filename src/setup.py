#!/usr/bin/env python

from distutils.core import setup

setup(name='sumac',
      version='2.1',
      description='SUMAC: Supermatrix Constructor',
      long_description='SUMAC (Supermatrix Constructor) is a Python package to data-mine GenBank, construct phylogenetic supermatrices, and assess the phylogenetic decisiveness of a matrix given the pattern of missing sequence data.',
      url='https://github.com/wf8/sumac',
      author='Will Freyman',
      license='GPLv3',
      keywords='phylogenetics supermatrix genbank',
      author_email='freyman@berkeley.edu',
      packages=['sumac'],
      package_dir={'sumac': ''},
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 2.7',
      ],
      #install_requires=['biopython'],
     )
