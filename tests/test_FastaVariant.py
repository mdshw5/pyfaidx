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



@pytest.fixture
def multiallelic_input(tmp_path):
    pysam = pytest.importorskip('pysam')
    pytest.importorskip('vcf')
    reference = tmp_path / 'reference.fa'
    reference.write_text('>test\nACGTTAGCCTAA\n>other\nACGT\n')
    vcf = tmp_path / 'variants.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=test,length=12>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHOM\tHET\tALT1\tREF\n'
        'test\t1\t.\tA\tC,G\t.\tPASS\t.\tGT\t2/2\t0|2\t1|2\t0/0\n'
        'test\t12\t.\tA\tC,G\t.\tPASS\t.\tGT\t2|2\t2|0\t2|1\t0/0\n'
    )
    bgzf = str(vcf) + '.gz'
    pysam.tabix_compress(str(vcf), bgzf)
    pysam.tabix_index(bgzf, preset='vcf')
    return str(reference), bgzf


@pytest.mark.parametrize('sample,expected,sites', [
    ('HOM', 'G', (1, 12)), ('HET', 'G', (1, 12)),
    ('ALT1', 'C', (1, 12)), ('REF', 'A', ()),
])
def test_variant_uses_genotype_alt_allele(multiallelic_input, sample, expected, sites):
    reference, vcf = multiallelic_input
    with FastaVariant(reference, vcf, sample=sample, as_raw=True) as fasta:
        assert fasta['test'][:1] == expected
        assert fasta['test'][-1:] == expected
        assert fasta['test'].variant_sites == sites


def test_variant_sites_contig_absent_from_vcf(multiallelic_input):
    reference, vcf = multiallelic_input
    with FastaVariant(reference, vcf, sample='HOM', as_raw=True) as fasta:
        assert fasta['other'][:] == 'ACGT'
        assert fasta['other'].variant_sites == ()
