import os
import pytest
from pyfaidx import Faidx, Fasta

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


@pytest.fixture
def remove_index():
    genes = Fasta('data/genes.fasta')
    del genes  # Support feature introduced in #111
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

def test_keys(remove_index):
    genes = Fasta('data/genes.fasta', key_function=get_gene_name)
    expect = ['BARD1', 'FGFR2', 'KF435149.1', 'MDM4', 'NR_104158.1', 'NR_104160.1', 'NR_104205.1', 'NR_104206.1', 'NR_104208.1', 'NR_104214.1', 'NR_104215.1', 'XR_241053.1', 'XR_241054.1', 'XR_241055.1', 'XR_241064.1', 'XR_241065.1', 'XR_241078.1', 'XR_241079.1', 'XR_241080.1', 'XR_241081.1']
    result = sorted(genes.keys())
    assert result == expect

def test_key_function_by_dictionary_get_key(remove_index):
    genes = Fasta('data/genes.fasta', key_function=get_gene_name)
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = genes['MDM4'][100-1:150]
    assert str(result) == expect

def test_key_function_by_fetch(remove_index):
    faidx = Faidx('data/genes.fasta', key_function=get_gene_name)
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = faidx.fetch('MDM4',
                         100, 150)
    assert str(result) == expect

def test_duplicated_keys(remove_index):
    with pytest.raises(ValueError):
        genes = Fasta('data/genes.fasta', key_function=get_duplicated_gene_name)

def test_duplicated_keys_shortest(remove_index):
    genes = Fasta('data/genes.fasta', key_function=get_duplicated_gene_name, duplicate_action="shortest")
    expect = 4573
    result = len(genes["BARD1"])
    assert expect == result

def test_duplicated_keys_longest(remove_index):
    genes = Fasta('data/genes.fasta', key_function=get_duplicated_gene_name, duplicate_action="longest")
    expect = 5317
    result = len(genes["BARD1"])
    assert expect == result
