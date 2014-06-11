import os
from pyfaidx import *
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)

ACCESSION_TO_GENE_NAME_DICT = {
    'AB821309.1': 'FGFR2',
    'KF435150.1': 'MDM4',
    'NR_104216.1': 'BARD1',
    # The rest are deliberately omitted
    # KF435149.1, NR_104215.1, NR_104212.1, NM_001282545.1 ...
    }

ACCESSION_TO_DUPLICATED_GENE_NAME_DICT = {
    'AB821309.1': 'FGFR2',
    'KF435150.1': 'MDM4',
    'NR_104216.1': 'BARD1',
    'NR_104215.1': 'BARD1', # Duplicated gene names will trigger a warning
    # The rest are deliberately omitted
    # KF435149.1, NR_104212.1, NM_001282545.1 ...
    }

def get_gene_name(accession):
    '''Return the gene name if found in ACCESSION_TO_GENE_NAME_DICT else return the original accession.'''
    return ACCESSION_TO_GENE_NAME_DICT.get(accession, accession)

def get_duplicated_gene_name(accession):
    '''Return the gene name if found in ACCESSION_TO_GENE_NAME_DICT else return the original accession.'''
    return ACCESSION_TO_DUPLICATED_GENE_NAME_DICT.get(accession, accession)


class TestFeatureKeyFunction:
    def __init__(self):
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.faidx = Faidx(self.fasta, key_function=get_gene_name)
        self.genes = Fasta(self.fasta, key_function=get_gene_name)

    def test_keys(self):
        expect = ['BARD1', 'FGFR2', 'KF435149.1', 'MDM4', 'NM_000465.3', 'NM_001282543.1', 'NM_001282545.1', 'NM_001282548.1', 'NM_001282549.1', 'NR_104212.1', 'NR_104215.1', 'XM_005249642.1', 'XM_005249643.1', 'XM_005249644.1', 'XM_005249645.1', 'XM_005265507.1', 'XM_005265508.1', 'XR_241079.1', 'XR_241080.1', 'XR_241081.1']
        result = sorted(self.genes.keys())
        assert result == expect

    def test_key_function_by_dictionary_get_key(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = self.genes['MDM4'][100-1:150]
        assert str(result) == expect

    def test_key_function_by_fetch(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = self.faidx.fetch('MDM4',
                             100, 150)
        assert str(result) == expect

    @raises(ValueError)
    def test_duplicated_keys(self):
        genes = Fasta(self.fasta, key_function=get_duplicated_gene_name)
