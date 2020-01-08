#!/usr/bin/env python
import os.path

ftp_timeout = 300

def fetch_genes(filename, suffix=None):
    from Bio import Entrez
    Entrez.email = "mdshw5@gmail.com"

    id_list = ['563317589', '557361099', '557361097', '543583796',
               '543583795', '543583794', '543583788', '543583786',
               '543583785', '543583740', '543583738', '530384540',
               '530384538', '530384536', '530384534', '530373237',
               '530373235', '530364726', '530364725', '530364724']

    search_results = Entrez.read(Entrez.epost("nucleotide", id=",".join(id_list)))

    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    records = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", webenv=webenv, query_key=query_key)

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

def fetch_chr22(filename):
    import gzip
    import shutil
    from io import BytesIO
    import ftplib

    remote = 'human_b36_male.fa.gz'
    ftp = ftplib.FTP('ftp-trace.ncbi.nih.gov', timeout=ftp_timeout) 
    ftp.login()
    ftp.cwd("1000genomes/ftp/pilot_data/technical/reference/")
    compressed = BytesIO()
    ftp.retrbinary('RETR %s' % remote, compressed.write)
    ftp.quit()
    compressed.seek(0)
    with gzip.GzipFile(fileobj = compressed) as gz:
        with open(filename, 'w') as fasta:
            chr22 = False
            for line in gz:
                if line[0:3] == '>22':
                    fasta.write(line)
                    chr22 = True
                elif not chr22:
                    continue
                elif chr22 and line[0] == '>':
                    curl.kill()
                    break
                elif chr22:
                    fasta.write(line)
    compressed.close()

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
    import ftplib

    ftp = ftplib.FTP('ftp-trace.ncbi.nih.gov', timeout=ftp_timeout) 
    ftp.login()
    ftp.cwd("1000genomes/ftp/pilot_data/release/2010_07/exon/snps/")
    with open(filename, 'wb') as vcf:
        ftp.retrbinary('RETR CEU.exon.2010_03.genotypes.vcf.gz', vcf.write)
    with open(filename, 'wb') as tbi:
        ftp.retrbinary('RETR CEU.exon.2010_03.genotypes.vcf.gz.tbi', tbi.write)
    ftp.quit()


if __name__ == "__main__":
    path = os.path.dirname(__file__)
    os.chdir(path)
    if not os.path.isfile("genes.fasta") or not os.path.isfile("genes.fasta.lower"):
        fetch_genes("genes.fasta")
    if not os.path.isfile("chr22.vcf.gz"):
        fetch_chr22_vcf("chr22.vcf.gz")
    if not os.path.isfile("chr22.fasta"):
        fetch_chr22("chr22.fasta")
    bgzip_compress_fasta("genes.fasta")
    bgzip_compress_fasta("chr22.fasta")
