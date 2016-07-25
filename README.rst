|Travis| |PyPI| |Landscape| |Coverage| |Depsy| |Appveyor|

Description
-----------

Samtools provides a function "faidx" (FAsta InDeX), which creates a
small flat index file ".fai" allowing for fast random access to any
subsequence in the indexed FASTA file, while loading a minimal amount of the
file in to memory. This python module implements pure Python classes for
indexing, retrieval, and in-place modification of FASTA files using a samtools
compatible index. The pyfaidx module is API compatible with the `pygr`_ seqdb module.
A command-line script "`faidx`_" is installed alongside the pyfaidx module, and
facilitates complex manipulation of FASTA files without any programming knowledge.

.. _`pygr`: https://github.com/cjlee112/pygr

If you use pyfaidx in your publication, please cite:

`Shirley MD`_, `Ma Z`_, `Pedersen B`_, `Wheelan S`_. `Efficient "pythonic" access to FASTA files using pyfaidx <https://dx.doi.org/10.7287/peerj.preprints.970v1>`_. PeerJ PrePrints 3:e1196. 2015.

.. _`Shirley MD`: http://github.com/mdshw5
.. _`Ma Z`: http://github.com/azalea
.. _`Pedersen B`: http://github.com/brentp
.. _`Wheelan S`: http://github.com/swheelan

Installation
------------

This package is tested under Linux, MacOS, and Windows using Python 3.2-3.4, 2.7, 2.6, and pypy and is available from the PyPI:

::

    pip install pyfaidx  # add --user if you don't have root

or download a `release <https://github.com/mdshw5/pyfaidx/releases>`_ and:

::

    python setup.py install

Usage
-----

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

    # Start attributes are 1-based
    >>> genes['NM_001282543.1'][200:230].start
    201

    # End attributes are 0-based
    >>> genes['NM_001282543.1'][200:230].end
    230

    >>> genes['NM_001282543.1'][200:230].longname
    'NM_001282543.1:201-230'

    >>> len(genes['NM_001282543.1'])
    5466

Note that start and end coordinates of Sequence objects are [1, 0]. This can be changed to [0, 0] by passing ``one_based_attributes=False`` to ``Fasta`` or ``Faidx``. This argument only affects the ``Sequence .start/.end`` attributes, and has no effect on slicing coordinates.

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

- Slicing start and end coordinates are 0-based, just like Python sequences.

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

Filter functions (returning True) limit the index:

.. code:: python

    # new in v0.3.8
    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta', filt_function = lambda x: x[0] == 'N')
    >>> genes.keys()
    dict_keys(['NR_104212', 'NM_001282543', 'NR_104216', 'NR_104215', 'NM_001282549', 'NM_000465', 'NM_001282545', 'NM_001282548'])
    >>> genes['XM_005249644']
    KeyError: XM_005249644 not in tests/data/genes.fasta.

Or just get a Python string:

.. code:: python

    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta', as_raw=True)
    >>> genes
    Fasta("tests/data/genes.fasta", as_raw=True)

    >>> genes['NM_001282543.1'][200:230]
    CTCGTTCCGCGCCCGCCATGGAACCGGATG

You can make sure that you always receive an uppercase sequence, even if your fasta file has lower case

.. code:: python

    >>> from pyfaidx import Fasta
    >>> reference = Fasta('tests/data/genes.fasta.lower', sequence_always_upper=True)
    >>> reference['gi|557361099|gb|KF435150.1|'][1:70]

    >gi|557361099|gb|KF435150.1|:2-70
    TGACATCATTTTCCACCTCTGCTCAGTGTTCAACATCTGACAGTGCTTGCAGGATCTCTCCTGGACAAA


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

Sequence names are truncated on any whitespace. This is a limitation of the indexing strategy. However, full names can be recovered:

.. code:: python

    # new in v0.3.7
    >>> from pyfaidx import Fasta
    >>> genes = Fasta('tests/data/genes.fasta')
    >>> for record in genes:
    ...   print(record.name)
    ...   print(record.long_name)
    ...
    gi|563317589|dbj|AB821309.1|
    gi|563317589|dbj|AB821309.1| Homo sapiens FGFR2-AHCYL1 mRNA for FGFR2-AHCYL1 fusion kinase protein, complete cds
    gi|557361099|gb|KF435150.1|
    gi|557361099|gb|KF435150.1| Homo sapiens MDM4 protein variant Y (MDM4) mRNA, complete cds, alternatively spliced
    gi|557361097|gb|KF435149.1|
    gi|557361097|gb|KF435149.1| Homo sapiens MDM4 protein variant G (MDM4) mRNA, complete cds
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

The FastaVariant class provides a way to integrate single nucleotide variant calls to generate a consensus sequence.

.. code:: python

    # new in v0.4.0
    >>> consensus = FastaVariant('tests/data/chr22.fasta', 'tests/data/chr22.vcf.gz', het=True, hom=True)
    RuntimeWarning: Using sample NA06984 genotypes.

    >>> consensus['22'].variant_sites
    (16042793, 21833121, 29153196, 29187373, 29187448, 29194610, 29821295, 29821332, 29993842, 32330460, 32352284)

    >>> consensus['22'][16042790:16042800]
    >22:16042791-16042800
    TCGTAGGACA

    >>> Fasta('tests/data/chr22.fasta')['22'][16042790:16042800]
    >22:16042791-16042800
    TCATAGGACA

    >>> consensus = FastaVariant('tests/data/chr22.fasta', 'tests/data/chr22.vcf.gz', sample='NA06984', het=True, hom=True, call_filter='GT == "0/1"')
    >>> consensus['22'].variant_sites
    (16042793, 29187373, 29187448, 29194610, 29821332)

.. _faidx:

It also provides a command-line script:

cli script: faidx
~~~~~~~~~~~~~~~~~

.. code:: bash

    Fetch sequences from FASTA. If no regions are specified, all entries in the
    input file are returned. Input FASTA file must be consistently line-wrapped,
    and line wrapping of output is based on input line lengths.

    positional arguments:
      fasta                 FASTA file
      regions               space separated regions of sequence to fetch e.g.
                            chr1:1-1000

    optional arguments:
      -h, --help            show this help message and exit
      -b BED, --bed BED     bed file of regions
      -o OUT, --out OUT     output file name (default: stdout)
      -i {bed,chromsizes,nucleotide,transposed}, --transform {bed,chromsizes,nucleotide,transposed} transform the requested regions into another format. default: None
      -c, --complement      complement the sequence. default: False
      -r, --reverse         reverse the sequence. default: False
      -a SIZE_RANGE, --size-range SIZE_RANGE
                            selected sequences are in the size range [low, high]. example: 1,1000 default: None
      -n, --no-names        omit sequence names from output. default: False
      -f, --full-names      output full names including description. default: False
      -x, --split-files     write each region to a separate file (names are derived from regions)
      -l, --lazy            fill in --default-seq for missing ranges. default: False
      -s DEFAULT_SEQ, --default-seq DEFAULT_SEQ
                            default base for missing positions and masking. default: N
      -d DELIMITER, --delimiter DELIMITER
                            delimiter for splitting names to multiple values (duplicate names will be discarded). default: None
      -g REGEX, --regex REGEX
                            selected sequences are those matching regular expression. default: .*
      -v, --invert-match    selected sequences are those not matching 'regions' argument. default: False
      -m, --mask-with-default-seq
                            mask the FASTA file using --default-seq default: False
      -M, --mask-by-case    mask the FASTA file by changing to lowercase. default: False
      --no-rebuild          do not rebuild the .fai index even if it is out of date. default: False
      --version             print pyfaidx version number

Examples:

.. code:: bash

    $ faidx tests/data/genes.fasta NM_001282543.1:201-210 NM_001282543.1:300-320
    >NM_001282543.1:201-210
    CTCGTTCCGC
    >NM_001282543.1:300-320
    GTAATTGTGTAAGTGACTGCA

    $ faidx --full-names tests/data/genes.fasta NM_001282543.1:201-210
    >NM_001282543.1| Homo sapiens BRCA1 associated RING domain 1 (BARD1), transcript variant 2, mRNA
    CTCGTTCCGC

    $ faidx --no-names tests/data/genes.fasta NM_001282543.1:201-210 NM_001282543.1:300-320
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

    $ faidx --regex "^NM_00128254[35]" genes.fasta
    >NM_001282543.1
    ..................
    ..................
    ..................
    >NM_001282545.1
    ..................
    ..................
    ..................

    $ faidx --lazy tests/data/genes.fasta NM_001282543.1:5460-5480
    >NM_001282543.1:5460-5480
    AAAAAAANNNNNNNNNNNNNN

    $ faidx --lazy --default-seq='Q' tests/data/genes.fasta NM_001282543.1:5460-5480
    >NM_001282543.1:5460-5480
    AAAAAAAQQQQQQQQQQQQQQ

    $ faidx tests/data/genes.fasta --bed regions.bed
    ...

    $ faidx --transform chromsizes tests/data/genes.fasta
    AB821309.1	3510
    KF435150.1	481
    KF435149.1	642
    NR_104216.1	4573
    NR_104215.1	5317
    NR_104212.1	5374
    ...

    $ faidx --transform bed tests/data/genes.fasta
    AB821309.1	1    3510
    KF435150.1	1    481
    KF435149.1	1    642
    NR_104216.1	1   4573
    NR_104215.1	1   5317
    NR_104212.1	1   5374
    ...

    $ faidx --transform nucleotide tests/data/genes.fasta
    name	start	end	A	T	C	G	N
    AB821309.1	1	3510	955	774	837	944	0
    KF435150.1	1	481	149	120	103	109	0
    KF435149.1	1	642	201	163	129	149	0
    NR_104216.1	1	4573	1294	1552	828	899	0
    NR_104215.1	1	5317	1567	1738	968	1044	0
    NR_104212.1	1	5374	1581	1756	977	1060	0
    ...

    faidx --transform transposed tests/data/genes.fasta
    AB821309.1	1	3510	ATGGTCAGCTGGGGTCGTTTCATC...
    KF435150.1	1	481	ATGACATCATTTTCCACCTCTGCT...
    KF435149.1	1	642	ATGACATCATTTTCCACCTCTGCT...
    NR_104216.1	1	4573	CCCCGCCCCTCTGGCGGCCCGCCG...
    NR_104215.1	1	5317	CCCCGCCCCTCTGGCGGCCCGCCG...
    NR_104212.1	1	5374	CCCCGCCCCTCTGGCGGCCCGCCG...
    ...

    $ faidx --split-files tests/data/genes.fasta
    $ ls
    AB821309.1.fasta	NM_001282549.1.fasta	XM_005249645.1.fasta
    KF435149.1.fasta	NR_104212.1.fasta	XM_005265507.1.fasta
    KF435150.1.fasta	NR_104215.1.fasta	XM_005265508.1.fasta
    NM_000465.3.fasta	NR_104216.1.fasta	XR_241079.1.fasta
    NM_001282543.1.fasta	XM_005249642.1.fasta	XR_241080.1.fasta
    NM_001282545.1.fasta	XM_005249643.1.fasta	XR_241081.1.fasta
    NM_001282548.1.fasta	XM_005249644.1.fasta

    $ faidx --delimiter='_' tests/data/genes.fasta 000465.3
    >000465.3
    CCCCGCCCCTCTGGCGGCCCGCCGTCCCAGACGCGGGAAGAGCTTGGCCGGTTTCGAGTCGCTGGCCTGC
    AGCTTCCCTGTGGTTTCCCGAGGCTTCCTTGCTTCCCGCTCTGCGAGGAGCCTTTCATCCGAAGGCGGGA
    .......



    $ faidx --size-range 5500,6000 -i chromsizes tests/data/genes.fasta
    NM_000465.3	5523

    $ faidx -m --bed regions.bed tests/data/genes.fasta
    ### Modifies tests/data/genes.fasta by masking regions using --default-seq character ###

    $ faidx -M --bed regions.bed tests/data/genes.fasta
    ### Modifies tests/data/genes.fasta by masking regions using lowercase characters ###


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


Changelog
---------

Please see the `releases <https://github.com/mdshw5/pyfaidx/releases>`_ for a
comprehensive list of version changes.

Known issues
------------

I try to fix as many bugs as possible, but most of this work is supported by a single developer. Please check the `known issues <https://github.com/mdshw5/pyfaidx/issues?utf8=✓&q=is%3Aissue+is%3Aopen+label%3Aknown>`_ for bugs relevant to your work. Pull requests are welcome.


Contributing
------------

Create a new Pull Request with one feauture. If you add a new feature, please
create also the relevant test.

To get test running on your machine:
 - Create a new virtualenv and install the `dev-requirements.txt`.
 - Download the test data running:

      python tests/data/download_gene_fasta.py

 - Run the tests with

      nosetests --with-coverage --cover-package=pyfaidx

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

.. |Coverage| image:: https://codecov.io/gh/mdshw5/pyfaidx/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/mdshw5/pyfaidx

.. |Depsy| image:: http://depsy.org/api/package/pypi/pyfaidx/badge.svg
   :target: http://depsy.org/package/python/pyfaidx

.. |Appveyor| image:: https://ci.appveyor.com/api/projects/status/80ihlw30a003596w?svg=true
   :target: https://ci.appveyor.com/project/mdshw5/pyfaidx
