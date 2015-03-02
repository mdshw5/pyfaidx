import os
from os.path import getmtime
from pyfaidx import Faidx, FastaIndexingError
from nose.tools import raises
from unittest import TestCase
from tempfile import NamedTemporaryFile
import time

path = os.path.dirname(__file__)
os.chdir(path)


class TestIndexing(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except FileNotFoundError:
            pass  # some tests may delete this file

    def test_build(self):
        expect_index = ("gi|563317589|dbj|AB821309.1|	3510	114	70	71\n"
        "gi|557361099|gb|KF435150.1|	481	3789	70	71\n"
        "gi|557361097|gb|KF435149.1|	642	4368	70	71\n"
        "gi|543583796|ref|NR_104216.1|	4573	5141	70	71\n"
        "gi|543583795|ref|NR_104215.1|	5317	9901	70	71\n"
        "gi|543583794|ref|NR_104212.1|	5374	15415	70	71\n"
        "gi|543583788|ref|NM_001282545.1|	4170	20980	70	71\n"
        "gi|543583786|ref|NM_001282543.1|	5466	25324	70	71\n"
        "gi|543583785|ref|NM_000465.3|	5523	30980	70	71\n"
        "gi|543583740|ref|NM_001282549.1|	3984	36696	70	71\n"
        "gi|543583738|ref|NM_001282548.1|	4113	40851	70	71\n"
        "gi|530384540|ref|XM_005249645.1|	2752	45151	70	71\n"
        "gi|530384538|ref|XM_005249644.1|	3004	48071	70	71\n"
        "gi|530384536|ref|XM_005249643.1|	3109	51246	70	71\n"
        "gi|530384534|ref|XM_005249642.1|	3097	54528	70	71\n"
        "gi|530373237|ref|XM_005265508.1|	2794	57830	70	71\n"
        "gi|530373235|ref|XM_005265507.1|	2848	60824	70	71\n"
        "gi|530364726|ref|XR_241081.1|	1009	63849	70	71\n"
        "gi|530364725|ref|XR_241080.1|	4884	65009	70	71\n"
        "gi|530364724|ref|XR_241079.1|	2819	70099	70	71\n")
        index_file = Faidx('data/genes.fasta').indexname
        result_index = open(index_file).read()
        assert result_index == expect_index

    def test_order(self):
        order = ("gi|563317589|dbj|AB821309.1|",
                 "gi|557361099|gb|KF435150.1|",
                 "gi|557361097|gb|KF435149.1|",
                 "gi|543583796|ref|NR_104216.1|",
                 "gi|543583795|ref|NR_104215.1|",
                 "gi|543583794|ref|NR_104212.1|",
                 "gi|543583788|ref|NM_001282545.1|",
                 "gi|543583786|ref|NM_001282543.1|",
                 "gi|543583785|ref|NM_000465.3|",
                 "gi|543583740|ref|NM_001282549.1|",
                 "gi|543583738|ref|NM_001282548.1|",
                 "gi|530384540|ref|XM_005249645.1|",
                 "gi|530384538|ref|XM_005249644.1|",
                 "gi|530384536|ref|XM_005249643.1|",
                 "gi|530384534|ref|XM_005249642.1|",
                 "gi|530373237|ref|XM_005265508.1|",
                 "gi|530373235|ref|XM_005265507.1|",
                 "gi|530364726|ref|XR_241081.1|",
                 "gi|530364725|ref|XR_241080.1|",
                 "gi|530364724|ref|XR_241079.1|")
        result = tuple(Faidx('data/genes.fasta').index.keys())
        assert result == order

    def test_valgrind_short_lines(self):
        """ Makes all full-length lines short and checks that error is raised
        in all appropriate circumstances.
        """
        indexed = []
        with open('data/genes.fasta') as genes:
            fasta = genes.readlines()
        n_lines = sum(1 for line in fasta)
        for n in range(n_lines):
            with NamedTemporaryFile(mode='w') as lines:
                for i, line in enumerate(fasta):
                    if i == n and line[0] != '>' and len(line) == 71:
                        line = line[:-3] + '\n'
                        full_line = True
                    elif i == n:
                        full_line = False
                    lines.write(line)
                    lines.flush()
                name = lines.name
                if full_line:
                    try:
                        Faidx(name)
                        indexed.append(True)
                    except FastaIndexingError:
                        indexed.append(False)
        assert not any(indexed)

    def test_valgrind_long_lines(self):
        """ Makes all full-length lines long and checks that error is raised
        in all appropriate circumstances.
        """
        indexed = []
        with open('data/genes.fasta') as genes:
            fasta = genes.readlines()
        n_lines = sum(1 for line in fasta)
        for n in range(n_lines):
            with NamedTemporaryFile(mode='w') as lines:
                for i, line in enumerate(fasta):
                    if i == n and line[0] != '>' and len(line) == 71:
                        line = line.rstrip('\n') + 'NNN' + '\n'
                        full_line = True
                    elif i == n:
                        full_line = False
                    lines.write(line)
                    lines.flush()
                name = lines.name
                if full_line:
                    try:
                        Faidx(name)
                        indexed.append(True)
                    except FastaIndexingError:
                        indexed.append(False)
        assert not any(indexed)

    def test_valgrind_blank_lines(self):
        """ Makes all full-length lines blank and checks that error is raised
        in all appropriate circumstances.
        """
        indexed = []
        with open('data/genes.fasta') as genes:
            fasta = genes.readlines()
        n_lines = sum(1 for line in fasta)
        for n in range(n_lines):
            with NamedTemporaryFile(mode='w') as lines:
                for i, line in enumerate(fasta):
                    if i == n and line[0] != '>' and len(line) == 71:
                        line = '\n'
                        full_line = True
                    elif i == n:
                        full_line = False
                    lines.write(line)
                    lines.flush()
                name = lines.name
                if full_line:
                    try:
                        Faidx(name)
                        indexed.append(True)
                    except FastaIndexingError:
                        indexed.append(False)
        assert not any(indexed)

    def test_reindex_on_modification(self):
        """ This test ensures that the index is regenerated when the FASTA
        modification time is newer than the index modification time.
        mdshw5/pyfaidx#50 """
        faidx = Faidx('data/genes.fasta')
        index_mtime = getmtime(faidx.indexname)
        faidx.close()
        os.utime('data/genes.fasta', (index_mtime + 10, ) * 2)
        time.sleep(2)
        faidx = Faidx('data/genes.fasta')
        assert getmtime(faidx.indexname) > index_mtime
