"""
Taken from https://github.com/brentp/pyfasta/blob/452d1ce5406ed73c4149b6d201bc65e4aa8afc27/tests/bench.py
"""

from itertools import islice
from tempfile import NamedTemporaryFile
import pyfaidx
import pyfasta
import pysam
from Bio import SeqIO
import time
import random
import os
import sys
from subprocess import call, check_output
import tracemalloc

random.seed(1234)

SEQLEN = 100000


def mean(s):
    return sum(s) / len(s)


def make_intervals(nreads=int(sys.argv[1]), seqlen=SEQLEN, readlen=int(sys.argv[2])):
    for _ in range(nreads):
        start = random.randint(0, seqlen)
        end = min(seqlen, start + readlen)
        yield (start, end)

intervals = tuple(make_intervals())


def make_long_fasta(filename, nrecs=2500, seqlen=SEQLEN):
    headers = []
    with open(filename, 'w') as f:
        s = "ACTGACTGAC"
        for i in range(nrecs):
            h = "header%i" % i
            headers.append(h)
            f.write('>' + h + '\n')
            for line in pyfaidx.wrap_sequence(80, s * (seqlen//10)):
                f.write(line)
    return headers

def read_dict(f, headers):
    for k in islice(headers, 10):
        for start, end in intervals:
            str(f[k][start:end])


def read_faidx(f, headers):
    for k in islice(headers, 10):
        for start, end in intervals:
            str(f.fetch(k, start + 1, end))


def read_fastahack(f, headers):
    for k in islice(headers, 10):
        for start, end in intervals:
            str(f.get_sub_sequence(k, start, end))


def read_pysam(f, headers):
    for k in islice(headers, 100):
        for start, end in intervals:
            str(pysam.faidx(f, '{0}:{1}-{2}'.format(k, start + 1, end)))


def read_samtools(f, headers):
    for k in islice(headers, 100):
        for start, end in intervals:
            check_output(['samtools', 'faidx', f, '{0}:{1}-{2}'.format(k, start + 1, end)])


def main():
    fa_file = NamedTemporaryFile()
    index = fa_file.name + '.fai'
    headers = make_long_fasta(fa_file.name)

    def pyfaidx_fasta(n):
        print('timings for pyfaidx.Fasta')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            f = pyfaidx.Fasta(fa_file.name)
            ti.append(time.time() - t)

            t = time.time()
            read_dict(f, headers)
            tf.append(time.time() - t)
            os.remove(index)
        # profile memory usage and report timings
        tracemalloc.start()
        f = pyfaidx.Fasta(fa_file.name)
        read_dict(f, headers)
        os.remove(index)
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    def pyfaidx_faidx(n):
        print('timings for pyfaidx.Faidx')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            f = pyfaidx.Faidx(fa_file.name)
            ti.append(time.time() - t)

            t = time.time()
            read_faidx(f, headers)
            tf.append(time.time() - t)
            os.remove(index)
        # profile memory usage and report timings
        tracemalloc.start()
        f = pyfaidx.Faidx(fa_file.name)
        read_faidx(f, headers)
        os.remove(index)
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    def fastahack_fetch(n):
        print('timings for fastahack.FastaHack')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            f = fastahack.FastaHack(fa_file.name)
            ti.append(time.time() - t)

            t = time.time()
            read_fastahack(f, headers)
            tf.append(time.time() - t)
            os.remove(index)
        # profile memory usage and report timings
        tracemalloc.start()
        f = fastahack.FastaHack(fa_file.name)
        read_fastahack(f, headers)
        os.remove(index)
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    def pyfasta_fseek(n):
        print('timings for pyfasta.Fasta (fseek)')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            f = pyfasta.Fasta(fa_file.name, record_class=pyfasta.FastaRecord)
            ti.append(time.time() - t)

            t = time.time()
            read_dict(f, headers)
            tf.append(time.time() - t)
            os.remove(fa_file.name + '.flat')
            os.remove(fa_file.name + '.gdx')
        # profile memory usage and report timings
        tracemalloc.start()
        f = pyfasta.Fasta(fa_file.name, record_class=pyfasta.FastaRecord)
        read_dict(f, headers)
        os.remove(fa_file.name + '.flat')
        os.remove(fa_file.name + '.gdx')
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    def pyfasta_fasta(n):
        print('timings for pyfasta.Fasta')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            f = pyfasta.Fasta(fa_file.name)
            ti.append(time.time() - t)

            t = time.time()
            read_dict(f, headers)
            tf.append(time.time() - t)
            os.remove(fa_file.name + '.flat')
            os.remove(fa_file.name + '.gdx')
        # profile memory usage and report timings
        tracemalloc.start()
        f = pyfasta.Fasta(fa_file.name)
        read_dict(f, headers)
        os.remove(fa_file.name + '.flat')
        os.remove(fa_file.name + '.gdx')
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    def pysam_faidx(n):
        print('timings for pysam.faidx')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            pysam.faidx(fa_file.name)
            ti.append(time.time() - t)

            t = time.time()
            read_pysam(fa_file.name, headers)
            tf.append(time.time() - t)
            os.remove(index)
        # profile memory usage and report timings
        tracemalloc.start()
        pysam.faidx(fa_file.name)
        read_pysam(fa_file.name, headers)
        os.remove(index)
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    def samtools_faidx(n):
        print('timings for samtools faidx')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            call(['samtools', 'faidx', fa_file.name])
            ti.append(time.time() - t)

            t = time.time()
            read_samtools(fa_file.name, headers)
            tf.append(time.time() - t)
            os.remove(index)
        print(mean(ti))
        print(mean(tf))

    def seqio_read(n):
        print('timings for Bio.SeqIO')
        ti = []
        tf = []
        for _ in range(n):
            t = time.time()
            fh = open(fa_file.name)
            f = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
            ti.append(time.time() - t)

            t = time.time()
            read_dict(f, headers)
            tf.append(time.time() - t)
            fh.close()
        # profile memory usage and report timings
        tracemalloc.start()
        fh = open(fa_file.name)
        f = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        read_dict(f, headers)
        fh.close()
        print(tracemalloc.get_traced_memory())
        print(mean(ti))
        print(mean(tf))
        tracemalloc.stop()

    n = 3
    pyfaidx_fasta(n)
    pyfaidx_faidx(n)
    pyfasta_fasta(n)
    pyfasta_fseek(n)
    seqio_read(n)
    samtools_faidx(1)
    pysam_faidx(1)


if __name__ == "__main__":
    main()
