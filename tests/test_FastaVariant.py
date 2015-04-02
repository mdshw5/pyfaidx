import os
from pyfaidx import FastaVariant, Fasta
from itertools import chain
from unittest import TestCase
from subprocess import call
import shlex

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaVariant(TestCase):

    def tearDown(self):
        try:
            os.remove('data/chr22.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_fetch_variant(self):
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        assert fasta['22'][32330458:32330462] == 'CAGG'  # het
        assert fasta['22'][32352282:32352286] == 'CAGC'  # hom

    def test_fetch_hom_variant(self):
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=False, as_raw=True)
        assert fasta['22'][32330458:32330462] == 'CGGG'  # het
        assert fasta['22'][32352282:32352286] == 'CAGC'  # hom

    def test_fetch_het_variant(self):
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=False, het=True, as_raw=True)
        assert fasta['22'][32330458:32330462] == 'CAGG'  # het
        assert fasta['22'][32352282:32352286] == 'CGGC'  # hom

    def test_all_pos(self):
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        assert fasta['22'].variant_sites == (16042793, 21833121, 29153196, 29187373, 29187448, 29194610, 29821295, 29821332, 29993842, 32330460, 32352284)

    def test_all_diff(self):
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        ref = Fasta('data/chr22.fasta', as_raw=True)
        assert all(ref['22'][pos-1] != fasta['22'][pos-1] for pos in fasta['22'].variant_sites)
