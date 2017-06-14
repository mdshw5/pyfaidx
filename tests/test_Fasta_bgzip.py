import os
from pyfaidx import Fasta, Faidx, UnsupportedCompressionFormat, FetchError
from itertools import chain
from unittest import TestCase
from nose.tools import raises
from nose.plugins.skip import SkipTest

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaBGZF(TestCase):
    def setUp(self):
        try:
            from Bio import SeqIO
        except ImportError:
            raise SkipTest

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.gz.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_integer_slice(self):
        fasta = Fasta('data/genes.fasta.gz')
        expect = fasta['gi|563317589|dbj|AB821309.1|'][:100].seq
        result = fasta[0][:100].seq
        assert expect == result

    def test_integer_index(self):
        fasta = Fasta('data/genes.fasta.gz')
        expect = fasta['gi|563317589|dbj|AB821309.1|'][100].seq
        result = fasta[0][100].seq
        assert expect == result

    def test_fetch_whole_fasta(self):
        expect = [line.rstrip('\n') for line in open('data/genes.fasta') if line[0] != '>']
        result = list(chain(*([line for line in record] for record in Fasta('data/genes.fasta.gz', as_raw=True))))
        assert expect == result

    def test_line_len(self):
        fasta = Fasta('data/genes.fasta.gz')
        for record in fasta:
            assert len(next(iter(record))) == fasta.faidx.index[record.name].lenc

    @raises(UnsupportedCompressionFormat)
    def test_mutable_bgzf(self):
        fasta = Fasta('data/genes.fasta.gz', mutable=True)

    def test_long_names(self):
        """ Test that deflines extracted using FastaRecord.long_name are
        identical to deflines in the actual file.
        """
        deflines = []
        with open('data/genes.fasta') as fasta_file:
            for line in fasta_file:
                if line[0] == '>':
                    deflines.append(line[1:-1])
        fasta = Fasta('data/genes.fasta.gz')
        long_names = []
        for record in fasta:
            long_names.append(record.long_name)
        assert deflines == long_names

    def test_fetch_whole_entry(self):
        faidx = Faidx('data/genes.fasta.gz')
        expect = ('ATGACATCATTTTCCACCTCTGCTCAGTGTTCAACATCTGA'
                'CAGTGCTTGCAGGATCTCTCCTGGACAAATCAATCAGGTACGACCA'
                'AAACTGCCGCTTTTGAAGATTTTGCATGCAGCAGGTGCGCAAGG'
                'TGAAATGTTCACTGTTAAAGAGGTCATGCACTATTTAGGTCAGTACAT'
                'AATGGTGAAGCAACTTTATGATCAGCAGGAGCAGCATATGGTATATTG'
                'TGGTGGAGATCTTTTGGGAGAACTACTGGGACGTCAGAGCTTCTCCGTG'
                'AAAGACCCAAGCCCTCTCTATGATATGCTAAGAAAGAATCTTGTCACTTT'
                'AGCCACTGCTACTACAGCAAAGTGCAGAGGAAAGTTCCACTTCCAGAAAAA'
                'GAACTACAGAAGACGATATCCCCACACTGCCTACCTCAGAGCATAAATGCA'
                'TACATTCTAGAGAAGGTGATTGAAGTGGGAAAAAATGATGACCTGGAGGACTC')
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             1, 481)
        assert str(result) == expect

    def test_fetch_middle(self):
        faidx = Faidx('data/genes.fasta.gz')
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             100, 150)
        assert str(result) == expect

    def test_fetch_end(self):
        faidx = Faidx('data/genes.fasta.gz')
        expect = 'TC'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 481)
        assert str(result) == expect

    @raises(FetchError)
    def test_fetch_border(self):
        """ Fetch past the end of a gene entry """
        faidx = Faidx('data/genes.fasta.gz')
        expect = 'TC'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 500)
        print(result)
        assert str(result) == expect

    def test_rev(self):
        faidx = Faidx('data/genes.fasta.gz')
        expect = 'GA'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 481)
        assert str(-result) == expect, result

    @raises(FetchError)
    def test_fetch_past_bounds(self):
        """ Fetch past the end of a gene entry """
        faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                                         480, 5000)

    @raises(FetchError)
    def test_fetch_negative(self):
        """ Fetch starting with a negative coordinate """
        faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                                         -10, 10)

    @raises(FetchError)
    def test_fetch_reversed_coordinates(self):
        """ Fetch starting with a negative coordinate """
        faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                                         50, 10)

    @raises(FetchError)
    def test_fetch_keyerror(self):
        """ Fetch a key that does not exist """
        faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
        result = faidx.fetch('gi|joe|gb|KF435150.1|',
                                         1, 10)

    def test_blank_string(self):
        """ seq[0:0] should return a blank string mdshw5/pyfaidx#53 """
        fasta = Fasta('data/genes.fasta.gz', as_raw=True)
        assert fasta['gi|557361099|gb|KF435150.1|'][0:0] == ''

    def test_slice_from_beginning(self):
        fasta = Fasta('data/genes.fasta.gz', as_raw=True)
        assert fasta['gi|557361099|gb|KF435150.1|'][:4] == 'ATGA'

    def test_slice_from_end(self):
        fasta = Fasta('data/genes.fasta.gz', as_raw=True)
        assert fasta['gi|557361099|gb|KF435150.1|'][-4:] == 'ACTC'

    def test_issue_74_start(self):
        f0 = Fasta('data/genes.fasta.gz', one_based_attributes=False)
        f1 = Fasta('data/genes.fasta.gz', one_based_attributes=True)
        assert f0['gi|557361099|gb|KF435150.1|'][0:90].start == f1['gi|557361099|gb|KF435150.1|'][0:90].start - 1

    def test_issue_74_consistency(self):
        f0 = Fasta('data/genes.fasta.gz', one_based_attributes=False)
        f1 = Fasta('data/genes.fasta.gz', one_based_attributes=True)
        assert str(f0['gi|557361099|gb|KF435150.1|'][0:90]) == str(f1['gi|557361099|gb|KF435150.1|'][0:90])

    def test_issue_74_end_faidx(self):
        f0 = Faidx('data/genes.fasta.gz', one_based_attributes=False)
        f1 = Faidx('data/genes.fasta.gz', one_based_attributes=True)
        end0 = f0.fetch('gi|557361099|gb|KF435150.1|', 1, 90).end
        end1 = f1.fetch('gi|557361099|gb|KF435150.1|', 1, 90).end
        assert end0 == end1

    def test_issue_74_end_fasta(self):
        f0 = Fasta('data/genes.fasta.gz', one_based_attributes=False)
        f1 = Fasta('data/genes.fasta.gz', one_based_attributes=True)
        end0 = f0['gi|557361099|gb|KF435150.1|'][1:90].end
        end1 = f1['gi|557361099|gb|KF435150.1|'][1:90].end
        print((end0, end1))
        assert end0 == end1

    def test_issue_79_fix(self):
        f = Fasta('data/genes.fasta.gz')
        s = f['gi|557361099|gb|KF435150.1|'][100:105]
        print((s.start, s.end))
        assert (101, 105) == (s.start, s.end)

    def test_issue_79_fix_negate(self):
        f = Fasta('data/genes.fasta.gz')
        s = f['gi|557361099|gb|KF435150.1|'][100:105]
        s = -s
        print((s.start, s.end))
        assert (105, 101) == (s.start, s.end)

    def test_issue_79_fix_one_based_false(self):
        f = Fasta('data/genes.fasta.gz', one_based_attributes=False)
        s = f['gi|557361099|gb|KF435150.1|'][100:105]
        print((s.start, s.end))
        assert (100, 105) == (s.start, s.end)

    def test_issue_79_fix_one_based_false_negate(self):
        f = Fasta('data/genes.fasta.gz', one_based_attributes=False)
        s = f['gi|557361099|gb|KF435150.1|'][100:105]
        print(s.__dict__)
        s = -s
        print(s.__dict__)
        assert (105, 100) == (s.start, s.end)

    @raises(UnsupportedCompressionFormat)
    def test_fetch_border_padded(self):
        """ Fetch past the end of a gene entry """
        faidx = Faidx('data/genes.fasta.gz', default_seq='N')
        expect = 'TCNNNNNNNNNNNNNNNNNNN'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 500)
        print(result)
        assert str(result) == expect
