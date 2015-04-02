#!/usr/bin/env python
import os.path

def fetch_genes(filename):
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
            if len(line) > 1:  # skip lines with only \n
                fasta.write(line)

def fetch_chr22(filename):
    from subprocess import Popen, PIPE

    grch36 = 'ftp://ftp-trace.ncbi.nih.gov//1000genomes/ftp/pilot_data/technical/reference/human_b36_male.fa.gz'
    curl = Popen(['curl', '-s', grch36], stdout=PIPE)
    gz = Popen(['gzip', '-dcq'], stdin=curl.stdout, stdout=PIPE)
    with gz.stdout as remote:
        with open(filename, 'w') as fasta:
            chr22 = False
            for line in remote:
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

def fetch_chr22_vcf(filename):
    from subprocess import call
    call(['curl', '-s', 'ftp://ftp-trace.ncbi.nih.gov//1000genomes/ftp/pilot_data/release/2010_07/exon/snps/CEU.exon.2010_03.genotypes.vcf.gz',
          '-o', filename])
    call(['curl', '-s', 'ftp://ftp-trace.ncbi.nih.gov//1000genomes/ftp/pilot_data/release/2010_07/exon/snps/CEU.exon.2010_03.genotypes.vcf.gz.tbi',
          '-o', filename + '.tbi'])



if __name__ == "__main__":
    path = os.path.dirname(__file__)
    os.chdir(path)
    fetch_genes("genes.fasta")
    fetch_chr22_vcf("chr22.vcf.gz")
    fetch_chr22("chr22.fasta")
