from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from bactopia import __version__

setup(
    author="Robert A. Petit III",
    author_email='robbie.petit@gmail.com',
    description="A Python package for working with Bactopia.",
    entry_points={
        'console_scripts': [
            'bactopia-summary=bactopia.commands.summary:main',
            'bactopia-jsonify=bactopia.commands.jsonify:main'
        ],
    },
    keywords=[],
    name='bactopia',
    packages=find_packages(),
    python_requires='>=3.6',
    url='https://github.com/bactopia/bactopia',
    version=__version__,
    zip_safe=False
)
