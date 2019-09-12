#!/usr/bin/env python3

from setuptools import setup, find_packages


setup(
    name='bactopia',
    version='1.0.1',
    description='An all-around workflow for bacterial genomics.',
    packages=find_packages(),
    package_data={},
    author='Robert A. Petit III',
    author_email='robert.petit@emory.edu',
    url='http://github.com/bactopia/bactopia',
    install_requires=[
        'ariba >= 2.14.1',
        'beautifulsoup4 >= 4.8.0',
        'biopython >= 1.74',
        'executor >= 21.3',
        'lxml >= 4.4.1',
        'ncbi-genome-download >= 0.2.10',
        'requests >= 2.22.0',
        'urllib3 >= 1.25.3',
    ],
    classifiers=[
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3 :: Only',
      'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    license='MIT',
)


