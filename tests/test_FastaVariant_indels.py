import os
import pytest

from pyfaidx import FastaVariant


path = os.path.dirname(__file__)
os.chdir(path)


@pytest.fixture(scope="module")
def ensure_tiny_vcf_tbi():
    """Compress and tabix-index tests/data/tiny.vcf -> tiny.vcf.gz + .tbi.
    Always reindex to pick up any edits to the VCF in this test module.
    """
    try:
        import pysam
    except ImportError:
        pytest.skip("pysam not installed.")

    vcf_txt = os.path.join("data", "tiny.vcf")
    vcf_gz = vcf_txt + ".gz"

    # (Re)create compressed and indexed VCF to ensure it's fresh
    if os.path.exists(vcf_gz):
        os.remove(vcf_gz)

    pysam.tabix_compress(vcf_txt, vcf_gz, force=True)

    if os.path.exists(vcf_gz + ".tbi"):
        os.remove(vcf_gz + ".tbi")

    pysam.tabix_index(vcf_gz, preset="vcf", force=True)

    assert os.path.exists(vcf_gz)
    assert os.path.exists(vcf_gz + ".tbi")
    return vcf_gz


def _make_fasta_variant(vcf_gz, **kwargs):
    return FastaVariant(
        os.path.join("data", "tiny.fa"), vcf_gz, as_raw=True, hom=True, het=True, **kwargs
    )


@pytest.mark.parametrize(
    "sample,exp_chr1,exp_chr2,exp_chr3",
    [
        (
            "sample1",
            # chr1 complete 48bp
            "AGAACCCCGGTTGGTTTAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            # chr2 complete 48bp
            "ACGTACGTANACGTACGTACGTACGTACGTACGTACGTACGT",
            # chr3 complete 20bp
            "GCGTGTTTACGTACGTACG",
        ),
        (
            "sample2",
            # chr1 first 48bp (the complete consensus sequence has two extra Ts at the end)
            "AAAACCCCGGTTGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTT",
            # chr2 complete 48bp
            "ACGTACGTANACGTACGTACGTACGTACGTACGTACGTACGT",
            # chr3 complete 20bp
            "ACGTACGTACGTACGTACGT",
        ),
        (
            "sample3",
            # chr1 complete 48bp
            "AGAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
            # chr2 complete 48bp
            "ACGTACGTANACGTACGTACGTACGTACGTACGTACGTACGT",
            # chr3 complete 20bp
            "GCGTACGTACGTACGTACGT",
        ),
    ],
)
def test_combined_indels_window(ensure_tiny_vcf_tbi, sample, exp_chr1, exp_chr2, exp_chr3):
    """Window 1..48 (0-based [0:48]) includes all variants on chr1, chr2 and chr3 for each sample.

    For sample1, this is equivalent to running:
        bcftools consensus -f tiny.fa -H A tiny.vcf.gz -s sample1
    """
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi, sample=sample)
        s1 = fasta["chr1"][0:48]
        assert s1 == exp_chr1

        s2 = fasta["chr2"][0:48]
        assert s2 == exp_chr2

        s3 = fasta["chr3"][0:20]
        assert s3 == exp_chr3
    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")


def test_chr1_deletion_only_window(ensure_tiny_vcf_tbi):
    """Window 12..19 (0-based [11:19]) includes only the deletion at pos15 (TTA->T)."""
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
        s = fasta["chr1"][11:19]
        assert s == "GTTTAA"  # 8 bp -> 6 bp (net -2)

    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")


def test_chr1_insertion_only_window(ensure_tiny_vcf_tbi):
    """Window 9..12 (0-based [8:12]) includes only the insertion at pos10 (G->GTT).
    Length remains 4 due to truncation to request size; content reflects insertion.
    """
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)

        s = fasta["chr1"][8:12]
        assert s == "GGTT"  # was 'GGGG'

    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")


def test_chr3_overlapping_del_then_ins(ensure_tiny_vcf_tbi):
    """Deletion TAC->T at pos4..6 (anchored at pos4), and insertion T->TTT at pos8.
    These are adjacent events. Validate consensus within a window spanning them.
    """
    fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
    s = fasta["chr3"][3:12]  # 0-based slice [3:12] == bases 4..12
    assert isinstance(s, str)
    assert len(s) == 9
    # Basic sanity: result starts with the anchor 'T' at pos4
    assert s[0] == "T"


def test_chr3_deletion_at_chrom_end(ensure_tiny_vcf_tbi):
    """End-of-chromosome deletion using anchored form (POS=19, REF=GT -> ALT=G).
    Ensures we handle trailing-edge edits without coordinate mismatch.
    """
    fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
    s = fasta["chr3"][12:20]
    assert s == "ACGTACG"  # one base shorter than reference
    assert len(s) == 7


def test_default_sample_is_first_with_warning(ensure_tiny_vcf_tbi):
    """When multiple samples exist and none specified, first sample is used with a warning."""
    try:
        with pytest.warns(RuntimeWarning):
            fv_default = _make_fasta_variant(ensure_tiny_vcf_tbi)
        fv_sample1 = _make_fasta_variant(ensure_tiny_vcf_tbi, sample="sample1")

        # Sanity: default should pick the first sample (sample1)
        s_def = fv_default["chr1"][0:12]
        s_s1 = fv_sample1["chr1"][0:12]
        assert s_def == s_s1
    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")
