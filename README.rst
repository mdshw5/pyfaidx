|Travis| |PyPI| |Landscape|

Description
-----------

Samtools provides a function "faidx" (FAsta InDeX), which creates a
small flat index file ".fai" allowing for fast random access to any
subsequence in the indexed FASTA file, while loading a minimal amount of the
file in to memory. This python module implements pure Python classes for
indexing, retrieval, and in-place modification of FASTA files.

A `manuscript <https://www.authorea.com/users/3882/articles/14770/_show_article>`_ is
currently under preparation.

Installation
------------

This package is tested under Python 3.2-3.4, 2.7, 2.6, and pypy.

::

    pip install pyfaidx

    or

    python setup.py install

Usage
-----

.. _benchmark: http://www.biostars.org/p/93364/#93390

.. code:: python

    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta')
    >>> genes
    Fasta("tests/data/genes.fasta")  # set strict_bounds=True for bounds checking

Acts like a dictionary.

.. code:: python

    >>> genes.keys() ('AB821309.1', 'KF435150.1', 'KF435149.1', 'NR_104216.1', 'NR_104215.1', 'NR_104212.1', 'NM_001282545.1', 'NM_001282543.1', 'NM_000465.3', 'NM_001282549.1', 'NM_001282548.1', 'XM_005249645.1', 'XM_005249644.1', 'XM_005249643.1', 'XM_005249642.1', 'XM_005265508.1', 'XM_005265507.1', 'XR_241081.1', 'XR_241080.1', 'XR_241079.1')

    >>> genes['NM_001282543.1'][200:230]
    >NM_001282543.1:201-230
    CTCGTTCCGCGCCCGCCATGGAACCGGATG

    >>> genes['NM_001282543.1'][200:230].seq
    'CTCGTTCCGCGCCCGCCATGGAACCGGATG'

    >>> genes['NM_001282543.1'][200:230].name
    'NM_001282543.1'

    >>> genes['NM_001282543.1'][200:230].start
    201

    >>> genes['NM_001282543.1'][200:230].end
    230

    >>> genes['NM_001282543.1'][200:230].longname
    'NM_001282543.1:201-230'

    >>> len(genes['NM_001282543.1'])
    5466

Indexes like a list:

.. code:: python

    >>> genes[0][:50]
    >AB821309.1:1-50
    ATGGTCAGCTGGGGTCGTTTCATCTGCCTGGTCGTGGTCACCATGGCAAC

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

Sequence can be buffered in memory using a read-ahead buffer
for fast sequential access:

.. code:: python

    >>> from timeit import timeit
    >>> fetch = "genes['NM_001282543.1'][200:230]"
    >>> read_ahead = "import pyfaidx; genes = pyfaidx.Fasta('tests/data/genes.fasta', read_ahead=10000)"
    >>> no_read_ahead = "import pyfaidx; genes = pyfaidx.Fasta('tests/data/genes.fasta')"
    >>> string_slicing = "genes = {}; genes['NM_001282543.1'] = 'N'*10000"

    >>> timeit(fetch, no_read_ahead, number=10000)
    0.2204863309962093
    >>> timeit(fetch, read_ahead, number=10000)
    0.1121859749982832
    >>> timeit(fetch, string_slicing, number=10000)
    0.0033553699977346696

Read-ahead buffering can reduce runtime by 1/2 for sequential accesses to buffered regions.

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

You can also perform line-based iteration, receiving the sequence lines as they appear in the FASTA file:

.. code:: python

    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta')
    >>> for line in genes['NM_001282543.1']:
    ...   print(line)
    CCCCGCCCCTCTGGCGGCCCGCCGTCCCAGACGCGGGAAGAGCTTGGCCGGTTTCGAGTCGCTGGCCTGC
    AGCTTCCCTGTGGTTTCCCGAGGCTTCCTTGCTTCCCGCTCTGCGAGGAGCCTTTCATCCGAAGGCGGGA
    CGATGCCGGATAATCGGCAGCCGAGGAACCGGCAGCCGAGGATCCGCTCCGGGAACGAGCCTCGTTCCGC
    ...

.. role:: red

If you want to modify the contents of your FASTA file in-place, you can use the `mutable` argument.
Any portion of the FastaRecord can be replaced with an equivalent-length string.
:red:`Warning`: *This will change the contents of your file immediately and permanently:*

.. code:: python

    >>> genes = Fasta('tests/data/genes.fasta', mutable=True)
    >>> type(genes['NM_001282543.1'])
    <class 'pyfaidx.MutableFastaRecord'>

    >>> genes['NM_001282543.1'][:10]
    >NM_001282543.1:1-10
    CCCCGCCCCT
    >>> genes['NM_001282543.1'][:10] = 'NNNNNNNNNN'
    >>> genes['NM_001282543.1'][:15]
    >NM_001282543.1:1-15
    NNNNNNNNNNCTGGC


It also provides a command-line script:

cli script: faidx
~~~~~~~~~~~~~~~~~

For usage type ``faidx -h``.

.. code:: bash

    $ faidx tests/data/genes.fasta NM_001282543.1:201-210 NM_001282543.1:300-320
    >NM_001282543.1:201-210
    CTCGTTCCGC
    >NM_001282543.1:300-320
    GTAATTGTGTAAGTGACTGCA

    $ faidx --no_names tests/data/genes.fasta NM_001282543.1:201-210 NM_001282543.1:300-320
    CTCGTTCCGC
    GTAATTGTGTAAGTGACTGCA

    $ faidx --complement tests/data/genes.fasta NM_001282543.1:201-210
    >NM_001282543.1:201-210 (complement)
    GAGCAAGGCG

    $ faidx --reverse tests/data/genes.fasta NM_001282543.1:201-210
    >NM_001282543.1:210-201
    CGCCTTGCTC

    $ faidx --reverse --complement tests/data/genes.fasta NM_001282543.1:201-210
    >NM_001282543.1:210-201 (complement)
    GCGGAACGAG

    $ faidx tests/data/genes.fasta NM_001282543.1
    >NM_001282543.1:1-5466
    CCCCGCCCCT........
    ..................
    ..................
    ..................

    $ faidx --lazy tests/data/genes.fasta NM_001282543.1:5460-5480
    >NM_001282543.1:5460-5480
    AAAAAAANNNNNNNNNNNNNN

    $ faidx --lazy --default_seq='Q' tests/data/genes.fasta NM_001282543.1:5460-5480
    >NM_001282543.1:5460-5480
    AAAAAAAQQQQQQQQQQQQQQ

    $ faidx tests/data/genes.fasta --bed regions.bed
    ...

    $ faidx --stats tests/data/genes.fasta
    AB821309.1	3510
    KF435150.1	481
    KF435149.1	642
    NR_104216.1	4573
    NR_104215.1	5317
    NR_104212.1	5374
    NM_001282545.1	4170
    NM_001282543.1	5466
    NM_000465.3	5523
    NM_001282549.1	3984
    NM_001282548.1	4113
    XM_005249645.1	2752
    XM_005249644.1	3004
    XM_005249643.1	3109
    XM_005249642.1	3097
    XM_005265508.1	2794
    XM_005265507.1	2848
    XR_241081.1	1009
    XR_241080.1	4884
    XR_241079.1	2819

    $ faidx --split_files tests/data/genes.fasta
    $ ls
    AB821309.1.fasta	NM_001282549.1.fasta	XM_005249645.1.fasta
    KF435149.1.fasta	NR_104212.1.fasta	XM_005265507.1.fasta
    KF435150.1.fasta	NR_104215.1.fasta	XM_005265508.1.fasta
    NM_000465.3.fasta	NR_104216.1.fasta	XR_241079.1.fasta
    NM_001282543.1.fasta	XM_005249642.1.fasta	XR_241080.1.fasta
    NM_001282545.1.fasta	XM_005249643.1.fasta	XR_241081.1.fasta
    NM_001282548.1.fasta	XM_005249644.1.fasta

Similar syntax as ``samtools faidx``


A lower-level Faidx class is also available:

.. code:: python

    >>> from pyfaidx import Faidx
    >>> fa = Faidx('genes.fa')  # can return str with as_raw=True
    >>> fa.index
    OrderedDict([('AB821309.1', IndexRecord(rlen=3510, offset=12, lenc=70, lenb=71)), ('KF435150.1', IndexRecord(rlen=481, offset=3585, lenc=70, lenb=71)),... ])

    >>> fa.index['AB821309.1'].rlen
    3510

    fa.fetch('AB821309.1', 1, 10)  # these are 1-based genomic coordinates
    >AB821309.1:1-10
    ATGGTCAGCT


-  If the FASTA file is not indexed, when ``Faidx`` is initialized the
   ``build_index`` method will automatically run, and
   the index will be written to "filename.fa.fai" with ``write_fai()``.
   where "filename.fa" is the original FASTA file.
-  Start and end coordinates are 1-based.


Changes
-------
*New in version 0.3.3*:

- `--split_files` option writes each returned sequence to an individual file. Names are generated based on the sequence name and region coordinates.
- `--stats` option prints the name and sequence length for each entry, suitable for use as a UCSC-style [chrom.sizes](http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes) file.
- Sequence `longname` attribute allows access to "chr:start-end (complement)" formatted names


*New in version 0.3.2*:

- Fasta `__getitem__` no longer initializes new FastaRecord classes
- Faidx `read_ahead` attribute implementaion avoids unnecessary disk hits (`#34 <https://github.com/mdshw5/pyfaidx/issues/34>`_)

*New in version 0.3.1*:

- Fasta can now accept an integer index in addition to string keys.

*New in version 0.3.0*:

- FastaRecord now works as a line-based iterator (`#30 <https://github.com/mdshw5/pyfaidx/issues/30>`_)
- Added MutableFastaRecord class that allows same-length in-place replacement for FASTA (`#29 <https://github.com/mdshw5/pyfaidx/issues/29>`_)

*New in version 0.2.9*:

- Added read-ahead buffer for fast sequential sequence access (`#26 <https://github.com/mdshw5/pyfaidx/issues/26>`_)
- Fixed a condition where `as_raw` parameter was not respected (`#27 <https://github.com/mdshw5/pyfaidx/issues/27>`_)

*New in version 0.2.8*:

- Small internal refactoring

*New in version 0.2.7*:

- Faidx and Fasta `strict_bounds` bounds checking logic is more correct
- Fasta `default_seq` parameter now works
- CLI script `faidx` now takes a BED file for fetching regions from a fasta

*New in version 0.2.6*:

- Faidx no longer has `raw_index` attribute or `rebuild_index` method (reduce memory footprint)
- Faidx index memory usage decreased by 31-40%
- .fai creation is streaming, performance increase for very large indices
- Possible speed regression when performing many small queries using `Fasta` class

*New in version 0.2.5*:

- Fasta and Faidx can take `default_seq` in addition to `as_raw`, `key_function`,
  and `strict_bounds` parameters.
- Fixed issue `#20 <https://github.com/mdshw5/pyfaidx/issues/20>`_
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
Shirley <http://mattshirley.com>`_, and was completed under the
mentorship and financial support of Drs. `Sarah
Wheelan <http://sjwheelan.som.jhmi.edu>`_ and `Vasan
Yegnasubramanian <http://yegnalab.onc.jhmi.edu>`_ at the Sidney Kimmel
Comprehensive Cancer Center in the Department of Oncology.

.. |Travis| image:: https://travis-ci.org/mdshw5/pyfaidx.svg?branch=master
    :target: https://travis-ci.org/mdshw5/pyfaidx

.. |PyPI| image:: https://img.shields.io/pypi/v/pyfaidx.svg?branch=master
    :target: https://pypi.python.org/pypi/pyfaidx

.. |Landscape| image:: https://landscape.io/github/mdshw5/pyfaidx/master/landscape.svg
   :target: https://landscape.io/github/mdshw5/pyfaidx/master
   :alt: Code Health
