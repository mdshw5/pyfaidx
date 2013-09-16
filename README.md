pyfaidx
=======

Installation
------------
This module is tested under Python 2.7, and will not install under any other version. Python 3 support is planned.

To install the `pyfaidx` module and `pyfaidx-fetch.py` script, simply run `sudo python setup.py install`.

Usage
-----

A pure python implementation of samtools faidx FASTA indexing.

    usage: pyfaidx-fetch.py [-h] [-r REGION] [-n] fasta
    
    Fetch sequence from faidx-indexed FASTA
    
    positional arguments:
      fasta                 faidx indexed FASTA file
    
    optional arguments:
      -h, --help            show this help message and exit
      -r REGION, --region REGION
                            region of sequence to fetch e.g. chr1:1-1000
      -n, --name            print sequence names

class faidx
-----------

    >>> from pyfaidx import faidx
    >>> fa = faidx('T7.fa')
    >>> fa.build('T7.fa', 'T7.fa.fai')
    >>> fa.index
    {'EM_PHG:V01146': {'lenc': 60, 'lenb': 61, 'rlen': 39937, 'offset': 40571}, 'EM_PHG:GU071091': {'lenc': 60, 'lenb': 61, 'rlen': 39778, 'offset': 74}}
    >>> fa.fetch('EM_PHG:V01146', 1, 10)
    EM_PHG:V01146
    TCTCACAGTG
    >>> fa.fetch('EM_PHG:V01146', 100, 120)
    EM_PHG:V01146
    GGTTGGGGATGACCCTTGGGT
    
- If the FASTA file is not indexed, when `faidx` is initialized the `build` method will automatically run,
producing "filename.fa.fai" where "filename.fa" is the original FASTA file.
- Start and end coordinates are 1-based.

Acknowledgements
------------------
This project is freely licensed by the author, [Matthew Shirley](http://mattshirley.com), and was completed under the mentorship 
and financial support of Drs. [Sarah Wheelan](http://sjwheelan.som.jhmi.edu) and [Vasan Yegnasubramanian](http://yegnalab.onc.jhmi.edu) at 
the Sidney Kimmel Comprehensive Cancer Center in the Department of Oncology.
