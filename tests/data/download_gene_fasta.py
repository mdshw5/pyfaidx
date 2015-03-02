#!/usr/bin/env python
import os.path

def fetch_fasta(filename):
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


if __name__ == "__main__":
    path = os.path.dirname(__file__)
    os.chdir(path)
    fasta_name = "genes.fasta"
    fetch_fasta(fasta_name)
