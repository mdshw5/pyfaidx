import os
import pytest
from pyfaidx import Fasta

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    fais = [
            "data/gene.bed12.fasta.fai",
            "data/chr17.hg19.part.fa.fai"
            ]
    for fai in fais:
        try:    
            os.remove(fai)
        except EnvironmentError:
            pass  # some tests may delete this file
    
def test_split_seq(remove_index):
    """ Fetch sequence by blocks """
    fa = Fasta('data/chr17.hg19.part.fa')
    
    gene = Fasta("data/gene.bed12.fasta")
    expect = gene[list(gene.keys())[0]][:].seq
    
    bed = "data/gene.bed12"
    with open(bed) as fi:
        record = fi.readline().strip().split("\t")

    chrom = record[0]
    start = int(record[1])
    strand = record[5]

    # parse bed12 format
    starts = [int(x) for x in record[11].split(",")[:-1]] 
    sizes = [int(x) for x in record[10].split(",")[:-1]]
    starts = [start + x  for x in starts]
    ends = [start + size  for start,size in zip(starts, sizes)] 
    
    # bed half-open
    if strand == "-":
        starts = [start + 1 for start in starts]
    else: 
        ends = [end - 1 for end in ends]
    
    intervals = zip(starts, ends) 
    result = fa.get_spliced_seq(chrom, intervals, rc=True)
    print(result.seq)
    print("====")
    print(expect)

    assert result.seq == expect

def test_get_seq_rc(remove_index):
    """ Check get_seq with rc argument """
    fa = Fasta('data/chr17.hg19.part.fa')
    
    result = fa.get_seq("chr17", 11, 20, rc=False)
    expect = "CCCTGTTCCT"
    print("normal")
    print(result.seq)
    print(expect)
    assert result.seq == expect
    
    result = fa.get_seq("chr17", 11, 20, rc=True)
    expect = "AGGAACAGGG"
    assert result.seq == expect
    print("rc")
    print(result.seq)
    print(expect)
