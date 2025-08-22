import os
import pytest
from pyfaidx import FastaVariant, Fasta

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/chr22.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

def test_fetch_variant(remove_index):
    try:
        import pysam
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        assert fasta['22'][32330458:32330462] == 'CAGG'  # het
        assert fasta['22'][32352282:32352286] == 'CAGC'  # hom
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")

def test_fetch_hom_variant(remove_index):
    try:
        import pysam
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=False, as_raw=True)
        assert fasta['22'][32330458:32330462] == 'CGGG'  # het
        assert fasta['22'][32352282:32352286] == 'CAGC'  # hom
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")

def test_fetch_het_variant(remove_index):
    try:
        import pysam
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=False, het=True, as_raw=True)
        assert fasta['22'][32330458:32330462] == 'CAGG'  # het
        assert fasta['22'][32352282:32352286] == 'CGGC'  # hom
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")

def test_fetch_chr_not_in_vcf(remove_index):
    try:
        import pysam
        fasta = FastaVariant('data/chr22andfake.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        assert fasta['fake'][:10] == 'ATCG' # fake is not in vcf 
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")
    
def test_all_pos(remove_index):
    try:
        import pysam
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        assert fasta['22'].variant_sites == (16042793, 21833121, 29153196, 29187373, 29187448, 29194610, 29821295, 29821332, 29993842, 32330460, 32352284)
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")

def test_all_diff(remove_index):
    try:
        import pysam
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', hom=True, het=True, as_raw=True)
        ref = Fasta('data/chr22.fasta', as_raw=True)
        print([(ref['22'][pos-1], fasta['22'][pos-1]) for pos in fasta['22'].variant_sites])
        assert all(ref['22'][pos-1] != fasta['22'][pos-1] for pos in fasta['22'].variant_sites)
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")


def test_call_filter_parameter(remove_index):
    try:
        import pysam
        # Valid call_filter
        fasta = FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', call_filter='GQ > 30', hom=True, het=True, as_raw=True)
        assert hasattr(fasta, 'filter')
        assert 'GQ' in fasta.filter and '30' in fasta.filter
        # Invalid call_filter (should raise ValueError)
        with pytest.raises(ValueError):
            FastaVariant('data/chr22.fasta', 'data/chr22.vcf.gz', call_filter='invalidfilter', hom=True, het=True, as_raw=True)
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")

def test_invalid_vcf_filename(remove_index):
    try:
        import pysam
        with pytest.raises((IOError, FileNotFoundError, OSError)):
            FastaVariant('data/chr22.fasta', 'data/does_not_exist.vcf.gz', hom=True, het=True, as_raw=True)
    except (ImportError, IOError):
        pytest.skip("pysam not installed.")

