import os
from pyfaidx import Fasta
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureSplicedSeq(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        fais = [
                "data/gene.bed12.fasta.fai",
                "data/chr17.hg19.part.fa.fai"
                ]
        for fai in fais:
            try:    
                os.remove(fai)
            except EnvironmentError:
                pass  # some tests may delete this file

    def test_split_seq(self):
        """ Fetch sequence by blocks """
        fa = Fasta('data/chr17.hg19.part.fa')
        
        gene = Fasta("data/gene.bed12.fasta")
        expect = gene[gene.keys()[0]][:].seq
        
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
        
        assert result.seq == expect
