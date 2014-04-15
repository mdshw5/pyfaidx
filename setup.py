from setuptools import setup

setup(
    name='pyfaidx',
    provides='pyfaidx',
    version="0.1.4",
    author='Matthew Shirley',
    author_email='mdshw5@gmail.com',
    url='http://mattshirley.com',
    description='"samtools faidx" compatible FASTA indexing in pure python',
    license='MIT',
    packages=['pyfaidx'],
    install_requires=['six'],
    entry_points={'console_scripts': ['faidx = pyfaidx.cli:main']},
    classifiers=[
            "Development Status :: 4 - Beta",
            "License :: OSI Approved :: MIT License",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: Unix",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.3",
            "Programming Language :: Python :: 3.2",
            "Programming Language :: Python :: 2.6",
            "Programming Language :: Python :: Implementation :: PyPy",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
