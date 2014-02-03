[![Build Status](https://travis-ci.org/mdshw5/pyfaidx.png?branch=master)](https://travis-ci.org/mdshw5/pyfaidx)

## Description

Samtools provides a function "faidx" (**fa**sta **i**n**d**e**x**), which 
creates a small flat index file ".fai" allowing for fast random access to any 
subsequence in the indexed fasta.

Pyfaidx provides an interface for creating and using this index for fast 
random access of subsequences in a "pythonic" manner. For example:

### class Genome

```python
>>> from pyfaidx import Genome
>>> genome = Genome('T7.fa')
>>> genome['EM_PHG:V01146'][0:10]
EM_PHG:V01146:1-10
TCTCACAGTG
```
    
It also provides a command-line script:

### cli script: pyfaidx

```shell
$ pyfaidx /tmp/hg19.fa -r chr10:1000000-1000010
GGAGGGCTGCA

$ pyfaidx /tmp/hg19.fa -n -r chr10:1000000-1000010
chr10:1000000-1000010
GGAGGGCTGCA
```

A lower-level Faidx class is also exposed:

### class Faidx
```python
>>> from pyfaidx import Faidx
>>> fa = Faidx('T7.fa')
>>> fa.build('T7.fa', 'T7.fa.fai')
>>> fa.index
{'EM_PHG:V01146': {'lenc': 60, 'lenb': 61, 'rlen': 39937, 'offset': 40571}, 'EM_PHG:GU071091': {'lenc': 60, 'lenb': 61, 'rlen': 39778, 'offset': 74}}
>>> fa.fetch('EM_PHG:V01146', 1, 10)
EM_PHG:V01146
TCTCACAGTG
>>> x = fa.fetch('EM_PHG:V01146', 100, 120)
>>> x
EM_PHG:V01146
GGTTGGGGATGACCCTTGGGT
>>> x.name
EM_PHG:V01146
>>> x.seq
GGTTGGGGATGACCCTTGGGT
```
    
- If the FASTA file is not indexed, when `faidx` is initialized the `build` method will automatically run,
producing "filename.fa.fai" where "filename.fa" is the original FASTA file.
- Start and end coordinates are 1-based.

## Installation

This package is tested under Python 3.3, 3.2, 2.7, 2.6, and pypy.

```
pip install -r requirements.txt
python setup.py install
```

## CLI Usage

"samtools faidx" compatible FASTA indexing in pure python.

    usage: pyfaidx [-h] [-r REGION] [-n] fasta
    
    Fetch sequence from faidx-indexed FASTA
    
    positional arguments:
      fasta                 faidx indexed FASTA file
    
    optional arguments:
      -h, --help            show this help message and exit
      -r REGION, --region REGION
                            region of sequence to fetch e.g. chr1:1-1000
      -n, --name            print sequence names

Acknowledgements
------------------
This project is freely licensed by the author, [Matthew Shirley](http://mattshirley.com), and was completed under the mentorship 
and financial support of Drs. [Sarah Wheelan](http://sjwheelan.som.jhmi.edu) and [Vasan Yegnasubramanian](http://yegnalab.onc.jhmi.edu) at 
the Sidney Kimmel Comprehensive Cancer Center in the Department of Oncology. Genome and Chromosome object implementations are influenced by 
the [Counsyl HGVS module](https://github.com/counsyl/hgvs).
