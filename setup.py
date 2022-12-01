from setuptools import setup
from io import open

install_requires = ['six', 'setuptools >= 0.7']

setup(
    name='pyfaidx',
    provides='pyfaidx',
    author='Matthew Shirley',
    author_email='mdshw5@gmail.com',
    url='http://mattshirley.com',
    description='pyfaidx: efficient pythonic random '
                'access to fasta subsequences',
    long_description=open('README.rst', encoding='utf-8').read(),
    license='BSD',
    packages=['pyfaidx'],
    install_requires=install_requires,
    use_scm_version={"local_scheme": "no-local-version"},
    setup_requires=['setuptools_scm'],
    entry_points={'console_scripts': ['faidx = pyfaidx.cli:main']},
    classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: BSD License",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: Unix",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: Implementation :: PyPy",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
