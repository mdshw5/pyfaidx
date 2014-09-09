|Travis| |PyPI|

Please cite `Shirley, Matthew (2014): pyfaidx: efficient pythonic random
access to fasta subsequences. figshare. DOI:10.6084/m9.figshare.972933`__.

.. __: http://dx.doi.org/10.6084/m9.figshare.972933


Description
-----------

Samtools provides a function "faidx" (FAsta InDeX), which creates a
small flat index file ".fai" allowing for fast random access to any
subsequence in the indexed fasta, while loading a minimal amount of the
file in to memory.

Pyfaidx provides an interface for creating and using this index for fast
random access of **DNA** subsequences from huge fasta files in a
"pythonic" manner. Indexing speed is comparable to samtools, and in some
cases sequence retrieval is much faster (benchmark_). For example:

.. _benchmark: http://www.biostars.org/p/93364/#93390

.. code:: python

    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta')
    >>> genes
    Fasta("tests/data/genes.fasta")  # set strict_bounds=True for bounds checking

Acts like a dictionary.

.. code:: python

    >>> genes.keys() ['NR_104215.1',
    'KF435150.1', 'NM_001282548.1', 'NM_001282549.1', 'XM_005249644.1',
    'NM_001282543.1', 'NR_104216.1', 'XM_005265508.1', 'XR_241079.1',
    'AB821309.1', 'XM_005249645.1', 'XR_241081.1', 'XM_005249643.1',
    'XM_005249642.1', 'NM_001282545.1', 'NR_104212.1', 'XR_241080.1',
    'XM_005265507.1', 'KF435149.1', 'NM_000465.3']

    >>> genes['NM_001282543.1'][200:230]
    >NM_001282543.1:201-230
    CTCGTTCCGCGCCCGCCATGGAACCGGATG

    >>> genes['NM_001282543.1'][200:230].seq
    'CTCGTTCCGCGCCCGCCATGGAACCGGATG'

    >>> genes['NM_001282543.1'][200:230].name
    'NM_001282543.1:201-230'

    >>> genes['NM_001282543.1'][200:230].start
    201

    >>> genes['NM_001282543.1'][200:230].end
    230

    >>> len(genes['NM_001282543.1'])
    5466

Slices just like a string:

.. code:: python

    >>> genes['NM_001282543.1'][200:230][:10]
    >NM_001282543.1:201-210
    CTCGTTCCGC

    >>> genes['NM_001282543.1'][200:230][::-1]
    >NM_001282543.1:230-201
    GTAGGCCAAGGTACCGCCCGCGCCTTGCTC

    >>> genes['NM_001282543.1'][200:230][::3]
    >NM_001282543.1:201-230
    CGCCCCTACA

    >>> genes['NM_001282543.1'][:]
    >NM_001282543.1:1-5466
    CCCCGCCCCT........

- Start and end coordinates are 0-based, just like Python.

Complements and reverse complements just like DNA

.. code:: python

    >>> genes['NM_001282543.1'][200:230].complement
    >NM_001282543.1 (complement):201-230
    GAGCAAGGCGCGGGCGGTACCTTGGCCTAC

    >>> genes['NM_001282543.1'][200:230].reverse
    >NM_001282543.1:230-201
    GTAGGCCAAGGTACCGCCCGCGCCTTGCTC

    >>> -genes['NM_001282543.1'][200:230]
    >NM_001282543.1 (complement):230-201
    CATCCGGTTCCATGGCGGGCGCGGAACGAG

Custom key functions provide cleaner access:

.. code:: python

    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta', key_function = lambda x: x.split('.')[0])
    >>> genes.keys()
    dict_keys(['NR_104212', 'NM_001282543', 'XM_005249644', 'XM_005249645', 'NR_104216', 'XM_005249643', 'NR_104215', 'KF435150', 'AB821309', 'NM_001282549', 'XR_241081', 'KF435149', 'XR_241079', 'NM_000465', 'XM_005265508', 'XR_241080', 'XM_005249642', 'NM_001282545', 'XM_005265507', 'NM_001282548'])
    >>> genes['NR_104212'][:10]
    >NR_104212:1-10
    CCCCGCCCCT

Or just get a Python string:

.. code:: python

    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta', as_raw=True)
    >>> genes
    Fasta("tests/data/genes.fasta", as_raw=True)

    >>> genes['NM_001282543.1'][200:230]
    CTCGTTCCGCGCCCGCCATGGAACCGGATG

It also provides a command-line script:

cli script: faidx
~~~~~~~~~~~~~~~~~

.. code:: shell

    $ faidx tests/data/genes.fasta NM_001282543.1:201-210 NM_001282543.1:300-320
    >NM_001282543.1
    CTCGTTCCGC
    >NM_001282543.1
    GTAATTGTGTAAGTGACTGCA

    $ faidx --complement tests/data/genes.fasta NM_001282543.1:201-210
    >NM_001282543.1
    GAGCAAGGCG

    $ faidx --reverse tests/data/genes.fasta NM_001282543.1:201-210
    >NM_001282543.1
    CGCCTTGCTC

    $ faidx tests/data/genes.fasta NM_001282543.1
    >NM_001282543.1
    CCCCGCCCCT........

    $ faidx tests/data/genes.fasta --list regions.txt
    ...

Similar syntax as ``samtools faidx``


A lower-level Faidx class is also available:

.. code:: python

    >>> from pyfaidx import Faidx
    >>> fa = Faidx('genes.fa')  # can return str with as_raw=True
    >>> fa.index
    OrderedDict([('AB821309.1', IndexRecord(rlen=3510, offset=12, lenc=70, lenb=71)), ('KF435150.1', IndexRecord(rlen=481, offset=3585, lenc=70, lenb=71)),... ])

    >>> fa.index['AB821309.1'].rlen
    3510

    fa.fetch('AB821309.1', 1, 10)
    >AB821309.1:1-10
    ATGGTCAGCT


-  If the FASTA file is not indexed, when ``Faidx`` is initialized the
   ``build_index`` method will automatically run, and
   the index will be written to "filename.fa.fai" with ``write_fai()``.
   where "filename.fa" is the original FASTA file.
-  Start and end coordinates are 1-based.

Installation
------------

This package is tested under Python 3.4, 3.3, 2.7, 2.6, and pypy.

::

    pip install pyfaidx

    or

    python setup.py install

CLI Usage
---------

::

    usage: faidx [-h] [-b BED] [-n] [--default_seq DEFAULT_SEQ] [--lazy]
                 [--complement] [--reverse]
                 fasta [regions [regions ...]]

    Fetch sequence from faidx-indexed FASTA

    positional arguments:
      fasta                 FASTA file
      regions               space separated regions of sequence to fetch e.g.
                            chr1:1-1000

    optional arguments:
      -h, --help            show this help message and exit
      -b BED, --bed BED     bed file of regions
      -n, --name            print sequence names. default: True
      --default_seq DEFAULT_SEQ
                            default base for missing positions. default: N
      --lazy                lazy region bounds checking - fill in default_seq for
                            missing ranges. default: False
      --complement          comlement the sequence. default: False
      --reverse             reverse the sequence. default: False

Changes
-------

*New in version 0.2.7*:

- Faidx and Fasta `strict_bounds` bounds checking logic is more correct
- Fasta `default_seq` parameter now works
- CLI script `faidx` now takes a BED file for fetching regions from a fasta

*New in version 0.2.6*:

- Faidx no longer has `raw_index` attribute or `rebuild_index` method (reduce memory footprint)
- Faidx index memory usage decreased by 31-40%
- *.fai creation is streaming, performance increase for very large indices
- Possible speed regression when performing many small queries using `Fasta` class

*New in version 0.2.5*:

- Fasta and Faidx can take `default_seq` in addition to `as_raw`, `key_function`,
  and `strict_bounds` parameters.
- Fixed issue `#20 <https://github.com/mdshw5/pyfaidx/issues/20>`__
- Faidx has attribute `raw_index` which is a list representing the fai file.
- Faidx has `rebuild_index` and `write_fai` functions for building and writing
  `raw_index` to file.
- Extra test cases, and test cases against Biopython SeqIO

*New in version 0.2.4*:

- Faidx index order is stable and non-random

*New in version 0.2.3*:

- Fixed a bug affecting Python 2.6

*New in version 0.2.2*:

- `Fasta` can receive the `strict_bounds` argument

*New in version 0.2.1*:

- `FastaRecord` str attribute returns a string
- `Fasta` is now an iterator

*New in version 0.2.0*:

- `as_raw` keyword arg for `Faidx` and `Fasta` allows a simple string return type
- `__str__` method for `FastaRecord` returns entire contig sequence

*New in version 0.1.9*:

- line wrapping of ``faidx`` is set based on the wrapping of the indexed
  fasta file
- added ``--reverse`` and ``--complement`` arguments to ``faidx``

*New in version 0.1.8*:

- ``key_function`` keyword argument to ``Fasta`` allows lookup based on function
  output

Acknowledgements
----------------

This project is freely licensed by the author, `Matthew
Shirley <http://mattshirley.com>`__, and was completed under the
mentorship and financial support of Drs. `Sarah
Wheelan <http://sjwheelan.som.jhmi.edu>`__ and `Vasan
Yegnasubramanian <http://yegnalab.onc.jhmi.edu>`__ at the Sidney Kimmel
Comprehensive Cancer Center in the Department of Oncology.

.. |Travis| image:: https://travis-ci.org/mdshw5/pyfaidx.svg?branch=master
    :target: https://travis-ci.org/mdshw5/pyfaidx

.. |PyPI| image:: https://img.shields.io/pypi/v/pyfaidx.svg?branch=master
    :target: https://pypi.python.org/pypi/pyfaidx
