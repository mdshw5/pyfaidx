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
    try:
        import vcf
    except ImportError:
        pytest.skip("PyVCF3 not installed.")

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
    kwargs.setdefault("indels", True)
    return FastaVariant(
        os.path.join("data", "tiny.fa"), vcf_gz,
        as_raw=True, hom=True, het=True, **kwargs
    )


def test_indels_disabled_by_default(ensure_tiny_vcf_tbi):
    """Without indels=True only SNPs are applied, as in previous releases."""
    with FastaVariant(
        os.path.join("data", "tiny.fa"), ensure_tiny_vcf_tbi,
        sample="sample1", as_raw=True,
    ) as fasta:
        # SNP at chr1:2 is applied; the chr1 indels are not.
        assert fasta["chr1"][:] == "AGAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"
        assert fasta["chr1"][8:12] == "GGGG"
        # The deletion at chr3:19 is not applied either.
        assert fasta["chr3"][12:20] == "ACGTACGT"


def test_mnp_not_applied_by_default(boundary_variants):
    with boundary_variants("2/2", as_raw=True, indels=False) as fasta:
        assert fasta.get_seq("chr", 3, 4) == "GT"


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
            # chr1 complete 50bp because the insertion at pos10 adds two bases
            "AAAACCCCGGTTGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
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
    The insertion expands the consensus without truncating the requested window.
    """
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)

        s = fasta["chr1"][8:12]
        assert s == "GGTTGG"

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


def test_chr3_deletion_overlapping_window_left_boundary(ensure_tiny_vcf_tbi):
    fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
    s = fasta["chr3"][2:5]  # 0-based slice [2:5] == bases 3..5
    assert s == "GT"


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


@pytest.fixture
def boundary_variants(tmp_path):
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("vcf")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr\nACGTACGTACGT\n>other\nACGT\n")
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr,length=12>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr\t4\t.\tTAC\tTA,GGC,TAACT\t.\tPASS\t.\tGT\t{gt}\n"
    )

    def make(gt, extra_records="", **kwargs):
        content = vcf.read_text().replace("{gt}", gt)
        content += extra_records
        path = tmp_path / ("variants_" + gt.replace("/", "_") + ".vcf")
        path.write_text(content)
        compressed = str(path) + ".gz"
        pysam.tabix_compress(str(path), compressed, force=True)
        pysam.tabix_index(compressed, preset="vcf", force=True)
        kwargs.setdefault("indels", True)
        return FastaVariant(str(reference), compressed, sample="S", **kwargs)

    return make


@pytest.mark.parametrize("gt,start,end,expected", [
    ("1/1", 5, 6, "A"),
    ("1/1", 6, 8, "GT"),
    ("1/1", 3, 5, "GTA"),
    ("2/2", 5, 6, "GC"),
    ("2/2", 3, 4, "GG"),
    ("2/1", 5, 6, "A"),
    ("3/3", 5, 6, "AACT"),
    ("3/3", 4, 5, "TA"),
    ("3/3", 6, 7, "ACTG"),
])
def test_edits_clipped_to_reference_window(boundary_variants, gt, start, end, expected):
    with boundary_variants(gt, as_raw=True) as fasta:
        assert fasta.get_seq("chr", start, end) == expected


def test_balanced_indels_clear_coordinates(ensure_tiny_vcf_tbi):
    with FastaVariant(os.path.join("data", "tiny.fa"), ensure_tiny_vcf_tbi,
                      sample="sample1", indels=True) as fasta:
        seq = fasta.get_seq("chr3", 4, 12)
        assert str(seq) == "TGTTTACGT"
        assert seq.start is None and seq.end is None
        assert seq[1].start is None


def test_overlapping_edit_is_skipped(boundary_variants):
    extra = "chr\t5\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n"
    with boundary_variants("1/1", extra_records=extra, as_raw=True) as fasta:
        assert fasta.get_seq("chr", 1, 8) == "ACGTAGT"
        assert fasta.get_seq("chr", 5, 8) == "AGT"


def test_mnp_keeps_coordinates_but_partial_indel_does_not(boundary_variants):
    with boundary_variants("2/2") as fasta:
        seq = fasta.get_seq("chr", 3, 4)
        assert str(seq) == "GG"
        assert (seq.start, seq.end) == (3, 4)
    with boundary_variants("1/1") as fasta:
        seq = fasta.get_seq("chr", 5, 5)
        assert str(seq) == "A"
        assert seq.start is None and seq.end is None


def test_clipped_reference_keeps_coordinates(boundary_variants):
    with boundary_variants("0/0") as fasta:
        seq = fasta.get_seq("chr", 10, 20)
        assert str(seq) == "CGT"
        assert (seq.start, seq.end) == (10, 12)
        assert str(fasta.get_seq("other", 1, 4)) == "ACGT"


def test_many_deletions_in_reference_order(tmp_path):
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("vcf")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr\n" + "ACGT" * 2000 + "\n")
    vcf = tmp_path / "variants.vcf"
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr,length=8000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    records = "".join(
        "chr\t{}\t.\tCG\tC\t.\tPASS\t.\tGT\t1/1\n".format(2 + 4 * i)
        for i in range(2000)
    )
    vcf.write_text(header + records)
    compressed = str(vcf) + ".gz"
    pysam.tabix_compress(str(vcf), compressed)
    pysam.tabix_index(compressed, preset="vcf")
    with FastaVariant(str(reference), compressed, sample="S", as_raw=True,
                      indels=True) as fasta:
        assert fasta["chr"][:] == "ACT" * 2000


@pytest.fixture
def overlap_variants(tmp_path):
    """Build a 40 bp chromosome and a two-record VCF on demand."""
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("vcf")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr\n" + "ACGT" * 10 + "\n")
    created = []
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr,length=40>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )

    def make(records, **kwargs):
        vcf = tmp_path / "variants_{}.vcf".format(len(created))
        created.append(vcf)
        vcf.write_text(
            header
            + "".join(
                "chr\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t1/1\n".format(pos, ref, alt)
                for pos, ref, alt in records
            )
        )
        compressed = str(vcf) + ".gz"
        pysam.tabix_compress(str(vcf), compressed, force=True)
        pysam.tabix_index(compressed, preset="vcf", force=True)
        kwargs.setdefault("indels", True)
        return FastaVariant(str(reference), compressed, sample="S", **kwargs)

    return make


def test_boundary_insertion_after_deletion(overlap_variants):
    """An insertion anchored on the final base of a deletion is applied.

    This matches ``bcftools consensus`` and is the pattern found in real
    chr22 1000 Genomes indel calls.
    """
    records = [(11, "GTACGTACGT", "G"), (20, "T", "TG")]
    with overlap_variants(records, as_raw=True) as fasta:
        # Whole chromosome, and windows wholly or partly inside the deletion.
        assert fasta["chr"][:] == "ACGTACGTACGGACGTACGTACGTACGTACGT"
        assert fasta["chr"][10:20] == "GG"
        assert fasta["chr"][19:20] == "G"
        assert fasta["chr"][11:19] == ""
        assert fasta["chr"][9:21] == "CGGA"


def test_boundary_insertion_extends_ref_span(overlap_variants):
    """A multi-base REF insertion at the boundary emits its remaining ref bases."""
    records = [(11, "GTACGTACGT", "G"), (20, "TAC", "TACG")]
    with overlap_variants(records, as_raw=True) as fasta:
        assert fasta["chr"][:] == "ACGTACGTACGACGGTACGTACGTACGTACGT"
        assert fasta["chr"][19:20] == ""
        assert fasta["chr"][9:21] == "CGA"
        assert fasta["chr"][19:22] == "ACG"


def test_boundary_deletion_after_deletion(overlap_variants):
    records = [(11, "GTACGTACGT", "G"), (20, "TAC", "T")]
    with overlap_variants(records, as_raw=True) as fasta:
        assert fasta["chr"][:] == "ACGTACGTACGGTACGTACGTACGTACGT"
        assert fasta["chr"][19:20] == ""
        assert fasta["chr"][20:25] == "GTA"
        assert fasta["chr"][9:21] == "CG"


def test_boundary_deletion_with_multibase_alt(overlap_variants):
    """A boundary deletion retains the ALT bases that map past the boundary."""
    records = [(11, "GTACGTACGT", "G"), (20, "TACG", "TA")]
    with overlap_variants(records, as_raw=True) as fasta:
        assert fasta["chr"][:] == "ACGTACGTACGATACGTACGTACGTACGT"
        assert fasta["chr"][19:20] == ""
        assert fasta["chr"][19:22] == "A"
        assert fasta["chr"][20:25] == "ATA"


def test_boundary_insertion_after_snp(overlap_variants):
    """A substitution at the boundary leaves its base in place for an insertion."""
    records = [(11, "G", "T"), (11, "G", "GG")]
    with overlap_variants(records, as_raw=True) as fasta:
        assert fasta["chr"][:] == "ACGTACGTACTGTACGTACGTACGTACGTACGTACGTACGT"
        assert fasta["chr"][10:20] == "TGTACGTACGT"


def test_boundary_indel_skipped_after_insertion(overlap_variants):
    """bcftools does not allow a boundary indel after a previous insertion."""
    records = [(11, "G", "GTT"), (11, "G", "GAA")]
    with overlap_variants(records, as_raw=True) as fasta:
        assert fasta["chr"][:] == "ACGTACGTACGTTTACGTACGTACGTACGTACGTACGTACGT"
        assert fasta["chr"][10:20] == "GTTTACGTACGT"


def test_non_anchored_overlap_is_skipped(overlap_variants):
    """A complex overlapping allele is still skipped, as in bcftools."""
    records = [(11, "GTACGTACGT", "G"), (20, "TA", "TGG")]
    with overlap_variants(records, as_raw=True) as fasta:
        assert fasta["chr"][:] == "ACGTACGTACGACGTACGTACGTACGTACGT"
        assert fasta["chr"][19:20] == ""
