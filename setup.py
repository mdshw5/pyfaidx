from distutils.core import setup
import sys
import os.path

__version__ = '0.0.1'

if sys.version_info < (2, 7):
        sys.exit("Python 2.7 is required.\n")

setup(
        name = 'pyfaidx',
        provides = 'pyfaidx',
        version = __version__,
        author = 'Matthew Shirley',
        author_email = 'mdshw5@gmail.com',
        url = 'http://mattshirley.com',
        description = 'A pure python implementation of samtools faidx FASTA indexing',
        license = 'MIT',
        py_modules = ['pyfaidx'],
        scripts = ['scripts/pyfaidx-fetch.py'],
        classifiers = [
                "Development Status :: 3 - Alpha",
                "License :: OSI Approved :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Research",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 2.7",
                "Topic :: Scientific/Engineering :: Bio-Informatics"
        ]
)
