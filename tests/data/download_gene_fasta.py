#!/usr/bin/env python
import os.path

def fetch_genes(filename, suffix=None):
    from Bio import Entrez
    Entrez.email = "mdshw5@gmail.com"

    # Use accession.version identifiers instead of GI numbers
    id_list = [
        'AB821309.1',  # was 563317589
        'KF435150.1',  # was 557361099
        'KF435149.1',  # was 557361097
        'NR_104216.1', # was 543583796
        'NR_104215.1', # was 543583795
        'NR_104214.1', # was 543583794
        'NR_104208.1', # was 543583788
        'NR_104206.1', # was 543583786
        'NR_104205.1', # was 543583785
        'NR_104160.1', # was 543583740
        'NR_104158.1', # was 543583738
        'XR_241081.1', # was 530384540
        'XR_241080.1', # was 530384538
        'XR_241079.1', # was 530384536
        'XR_241078.1', # was 530384534
        'XR_241065.1', # was 530373237
        'XR_241064.1', # was 530373235
        'XR_241055.1', # was 530364726
        'XR_241054.1', # was 530364725
        'XR_241053.1', # was 530364724
    ]

    # Directly fetch by accession.version (no EPost needed)
    records = Entrez.efetch(db="nucleotide", id=",".join(id_list), rettype="fasta", retmode="text")

    with open(filename, 'w') as fasta:
        for line in records:
            if suffix is not None and line[0] == '>':
                line = line.rstrip('\n')
                line = ''.join([line, suffix, '\n'])
            if len(line) == 1:  # skip lines with only \n
                continue
            fasta.write(line)
    with open(filename, 'r') as fasta:
        with open('.'.join([filename, 'lower']), 'w') as lower:
            for line in fasta:
                if line[0] != '>':
                    line = line.lower()
                lower.write(line)

def issue_141_genes(newlines, carriages):
    """ Creates a fasta file with the genes from issue #141 
    https://github.com/mdshw5/pyfaidx/issues/141"""
    # replace the newline characters with carriage returns
    # open newlines file and read lines
    with open(newlines, 'r') as newlines_file:
        lines = newlines_file.readlines()
    # open carriages file and write lines
    with open(carriages, 'w') as carriages_file:
        for line in lines:
            # replace newline with carriage return
            line = line.replace('\n', '\r\n')
            # write the modified line to the file
            carriages_file.write(line)

def fetch_chr22(filename):
    import requests
    import gzip
    import io

    with requests.get('https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/technical/reference/human_b36_male.fa.gz') as compressed:
        with open(filename, 'w') as fasta, gzip.GzipFile(fileobj=io.BytesIO(compressed.content)) as gz:
            chr22 = False
            for line in gz:
                line_content = line.decode()
                if line_content[0:3] == '>22':
                    fasta.write(line_content)
                    chr22 = True
                elif not chr22:
                    continue
                elif chr22 and line_content[0] == '>':
                    break
                elif chr22:
                    fasta.write(line_content)

def add_fake_chr(existing_fasta, filename):
    with open(filename, 'w') as fasta:
        with open(existing_fasta, 'r') as old_fa:
            for line in old_fa:
                fasta.write(line)
        fasta.write('>fake chromosome not in vcf\n')
        fasta.write('ATCG\n')

def fake_chr22(filename):
    """ Fake up some data """
    chr22_len = 49691432
    mod_70 = chr22_len % 70
    with open(filename, 'w') as fake_file:
        fake_file.write('>22 fake data for testing\n')
        while chr22_len > mod_70:
            fake_file.write('N' * 70 + '\n')
            chr22_len -= 70
        fake_file.write('N' * mod_70 + '\n')

def bgzip_compress_fasta(filename):
    from Bio.bgzf import BgzfWriter
    with BgzfWriter(filename=filename + '.gz') as compressed, open(filename, 'r') as fasta:
        for line in fasta:
            compressed.write(line)

def fetch_chr22_vcf(filename):
    import requests
    
    with requests.get('https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2010_07/exon/snps/CEU.exon.2010_03.genotypes.vcf.gz') as vcf:
        with open(filename, 'wb') as out:
            out.write(vcf.content)
    with requests.get('https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/release/2010_07/exon/snps/CEU.exon.2010_03.genotypes.vcf.gz.tbi') as tbi:
        with open(filename + '.tbi', 'wb') as out:
            out.write(tbi.content)


if __name__ == "__main__":
    path = os.path.dirname(__file__)
    os.chdir(path)
    if not os.path.isfile("genes.fasta") or not os.path.isfile("genes.fasta.lower"):
        print("GETTING genes")
        fetch_genes("genes.fasta")
    if not os.path.isfile("issue_141.fasta"):
        print("GETTING issue 141 genes")
        issue_141_genes("genes.fasta", "issue_141.fasta")
    if not os.path.isfile("chr22.vcf.gz"):
        print("GETTING vcf")
        fetch_chr22_vcf("chr22.vcf.gz")
    if not os.path.isfile("chr22.fasta"):
        print("GETTING chr22.fasta")
        fetch_chr22("chr22.fasta")
    if not os.path.isfile("chr22andfake.fasta"):
        print("adding fake chr")
        add_fake_chr("chr22.fasta", "chr22andfake.fasta")
    bgzip_compress_fasta("genes.fasta")
    bgzip_compress_fasta("chr22.fasta")