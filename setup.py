from setuptools import setup
from io import open
import sys

install_requires = ['six', 'setuptools >= 0.7']
if sys.version_info[0] == 2 and sys.version_info[1] == 6:
    install_requires.extend(['ordereddict', 'argparse'])


def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


setup(
    name='pyfaidx',
    provides='pyfaidx',
    version=get_version(open('pyfaidx/__init__.py', encoding='utf-8').read()),
    author='Matthew Shirley',
    author_email='mdshw5@gmail.com',
    url='http://mattshirley.com',
    description='pyfaidx: efficient pythonic random '
                'access to fasta subsequences',
    long_description=open('README.rst', encoding='utf-8').read(),
    license='BSD',
    packages=['pyfaidx'],
    install_requires=install_requires,
    entry_points={'console_scripts': ['faidx = pyfaidx.cli:main']},
    classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: BSD License",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: Unix",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.3",
            "Programming Language :: Python :: 3.2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 2.6",
            "Programming Language :: Python :: Implementation :: PyPy",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
