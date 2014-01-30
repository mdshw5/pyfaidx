from setuptools import setup, Extension
import distutils.core
import sys
import os
setup(
        name = 'pyfaidx',
        provides = 'pyfaidx',
        version = "0.1.0",
        author = 'Matthew Shirley',
        author_email = 'mdshw5@gmail.com',
        url = 'http://mattshirley.com',
        description = 'A pure python implementation of samtools faidx FASTA indexing',
        license = 'MIT',
        packages = ['pyfaidx'],
        entry_points = { 'console_scripts': [ 'pyfaidx = pyfaidx.cli:main' ] },
        classifiers = [
                "Development Status :: 3 - Alpha",
                "License :: OSI Approved :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Research",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 3.3",
                "Topic :: Scientific/Engineering :: Bio-Informatics"
        ]
)