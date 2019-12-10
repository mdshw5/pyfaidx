import os
import sys
from pyfaidx import Fasta
from tempfile import NamedTemporaryFile
from unittest import TestCase
from nose.tools import raises
from difflib import Differ

path = os.path.dirname(__file__)
os.chdir(path)


class TestFastaRecord(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_sequence_uppercase(self):
        """Test that the sequence is always returned in
        uppercase, even if it is in lowercase in the
        reference genome.
        """
        filename = "data/genes.fasta.lower"
        reference_upper = Fasta(filename, sequence_always_upper=True)
        reference_normal = Fasta(filename)
        os.remove('data/genes.fasta.lower.fai')
        assert reference_upper['gi|557361099|gb|KF435150.1|'][
            1:100].seq == reference_normal['gi|557361099|gb|KF435150.1|'][
                1:100].seq.upper()

    def test_long_names(self):
        """ Test that deflines extracted using FastaRecord.long_name are
        identical to deflines in the actual file.
        """
        deflines = []
        with open('data/genes.fasta') as fasta_file:
            for line in fasta_file:
                if line[0] == '>':
                    deflines.append(line[1:-1])
        fasta = Fasta('data/genes.fasta')
        long_names = []
        for record in fasta:
            long_names.append(record.long_name)
        print(tuple(zip(deflines, long_names)))
        assert deflines == long_names

    def test_issue_62(self):
        """ Check for pathogenic FastaRecord.long_name behavior in mdshw5/pyfaidx#62 """
        deflines = []
        line_len = None
        with open('data/genes.fasta', 'rb') as fasta_file:
            with open('data/issue_62.fa', 'wb') as fasta_uniform_len:
                for line in fasta_file:
                    if line.startswith(b'>'):
                        deflines.append(line[1:-1].decode('ascii'))
                        fasta_uniform_len.write(line)
                    elif line_len is None:
                        line_len = len(line)
                        fasta_uniform_len.write(line)
                    elif line_len > len(line):
                        fasta_uniform_len.write(line.rstrip() + b'N' *
                                                (line_len - len(line)) + b'\n')
                    else:
                        fasta_uniform_len.write(line)
        fasta = Fasta('data/issue_62.fa', as_raw=True)
        long_names = []
        for record in fasta:
            long_names.append(record.long_name)
        try:
            os.remove('data/issue_62.fa')
            os.remove('data/issue_62.fa.fai')
        except EnvironmentError:
            pass
        sys.stdout.writelines(tuple(Differ().compare(deflines, long_names)))
        assert deflines == long_names

    def test_unpadded_length(self):
        filename = "data/padded.fasta"
        with open(filename, 'w') as padded:
            padded.write(">test_padded\n")
            for n in range(10):
                padded.write("N" * 80)
                padded.write("\n")
            padded.write("N" * 30)
            padded.write("A" * 20)
            padded.write("N" * 30)
            padded.write("\n")
            for n in range(10):
                padded.write("N" * 80)
                padded.write("\n")

        fasta = Fasta(filename)
        expect = 20
        result = fasta["test_padded"].unpadded_len
        print(expect, result)
        assert expect == result
        os.remove('data/padded.fasta')
        os.remove('data/padded.fasta.fai')

    def test_numpy_array(self):
        """ Test the __array_interface__ """
        import numpy
        filename = "data/genes.fasta.lower"
        reference = Fasta(filename)
        np_array = numpy.asarray(reference[0])
        assert isinstance(np_array, numpy.ndarray)


class TestMutableFastaRecord(TestCase):
    def setUp(self):
        with open('data/genes_mutable.fasta', 'wb') as mutable:
            mutable.write(open('data/genes.fasta', 'rb').read())
        self.mutable_fasta = Fasta('data/genes_mutable.fasta', mutable=True)

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file
        try:
            os.remove('data/genes_mutable.fasta')
        except EnvironmentError:
            pass  # some tests may delete this file
        try:
            os.remove('data/genes_mutable.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_mutate_fasta_to_same(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=True)
        fasta = Fasta('data/genes.fasta', mutable=False)
        chunk = fasta['gi|557361099|gb|KF435150.1|'][0:100]
        mutable['gi|557361099|gb|KF435150.1|'][0:100] = chunk.seq
        assert str(fasta['gi|557361099|gb|KF435150.1|']) == str(
            mutable['gi|557361099|gb|KF435150.1|'])

    def test_mutate_fasta_to_N(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=True)
        chunk = 100 * 'N'
        mutable['gi|557361099|gb|KF435150.1|'][0:100] = chunk
        assert mutable['gi|557361099|gb|KF435150.1|'][0:100].seq == chunk

    def test_mutate_single_position(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=True)
        chunk = 'N'
        mutable['gi|557361099|gb|KF435150.1|'][0] = chunk
        assert mutable['gi|557361099|gb|KF435150.1|'][0].seq == chunk

    @raises(TypeError)
    def test_mutate_immutable_fasta(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=False)
        chunk = 100 * 'N'
        mutable['gi|557361099|gb|KF435150.1|'][0:100] = chunk

    @raises(IOError)
    def test_mutate_too_long(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=True)
        chunk = 101 * 'N'
        mutable['gi|557361099|gb|KF435150.1|'][0:100] = chunk
