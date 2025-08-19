import os
import pytest
from collections import OrderedDict
import threading
from pyfaidx import Fasta
import random
import tempfile
import time
import shutil
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)


class _ThreadReadSequence(threading.Thread):
    def __init__(self, rand, result_map, result_lock, name, seq):
        super(_ThreadReadSequence, self).__init__()

        seq_len = len(seq)
        sub_seq_slices = list(slice(i, min(i + 20, seq_len)) for i in range(0, seq_len, 20))
        random.shuffle(sub_seq_slices)

        self.result_map = result_map
        self.result_lock = result_lock
        self.name = name
        self.seq = seq
        self.sub_seq_slices = sub_seq_slices

    def run(self):
        name = self.name
        seq = self.seq

        sub_seqs = [''] * len(self.sub_seq_slices)
        for sub_seq_slice in self.sub_seq_slices:
            sub_seqs[sub_seq_slice.start//20] = seq[sub_seq_slice]
            time.sleep(0)

        # Put sub-sequences in correct order
        seq_str = ''.join(sub_seqs)

        with self.result_lock:
            self.result_map[name] = seq_str


class _ThreadWriteSequence(threading.Thread):
    def __init__(self, rand, name, seq):
        super(_ThreadWriteSequence, self).__init__()

        seq_len = len(seq)
        sub_seq_slices = list(slice(i, min(i + 20, seq_len)) for i in range(0, seq_len, 20))
        random.shuffle(sub_seq_slices)

        self.name = name
        self.seq = seq
        self.sub_seq_slices = sub_seq_slices

    def run(self):
        seq = self.seq
        seq_len = len(seq)
        seq_str = seq[:].lower()

        for sub_seq_slice in self.sub_seq_slices:
            try:
                seq[sub_seq_slice] = seq_str[sub_seq_slice]
                time.sleep(0)
            except Exception:
                # Conflicting simultaneous writes are likely to cause exceptions
                # We test for the expected string at the end, so ignore interim
                # failures.
                pass


class TestFastaIntIndex(TestCase):
    def setUp(self):
        self.longMessage = True
        self.maxDiff = None
        self.tmp_dir = tempfile.mkdtemp()
        # Use a seeded random orders are randomish within a test, but the same across test runs
        self.rand = random.Random(8903423147)

    def tearDown(self):
        tmp_dir = getattr(self, 'tmp_dir', None)
        if tmp_dir:
            shutil.rmtree(tmp_dir)

        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_simultaneous_reads(self):
        """
        Test that each read of a sequence range is atomic.
        To do this, spawn several threads to simultaneously read the sequences
        in a Fasta file in pieces. If the reads are not atomic, then it is
        reasonably likely (with sufficient concurrency) that a read from one
        thread will affect that in another, so the sequences will not be
        read properly.
        """
        # Read in original file data
        ref_seq_map = OrderedDict()
        with Fasta('data/genes.fasta', as_raw=True, strict_bounds=True) as fasta:
            for name, seq in fasta.records.items():
                ref_seq_map[name] = seq[:]

        # Initialize map with fasta sequence names to enforce same ordering as 'ref_seq_map'
        thread_result_lock = threading.Lock()
        thread_read_seq_map = OrderedDict((name, None) for name in ref_seq_map)

        # Read file again, using many threads and simultaneously reading each sequence in pieces
        with Fasta('data/genes.fasta', as_raw=True, strict_bounds=True) as fasta:
            threads = []
            for name, seq in fasta.records.items():
                threads.append(_ThreadReadSequence(self.rand, thread_read_seq_map, thread_result_lock, name, seq))

            for thread in threads:
                thread.start()

            for thread in threads:
                thread.join()

        self.assertEqual(thread_read_seq_map, ref_seq_map)

    def test_simultaneous_writes(self):
        """
        Test that each write of a sequence range is atomic.
        To do this, spawn several threads to simultaneously write sequences
        to a Fasta file in pieces. If the writes are not atomic, then it is
        reasonably likely (with sufficient concurrency) that a write from one
        thread will affect that in another, so the sequences will not be
        written properly. To make sure all sequences are mutated, the writes
        will transform the sequence to lower-case.
        """

        tmp_dir = self.tmp_dir

        tmp_fasta = os.path.join(tmp_dir, 'genes_write.fasta')
        shutil.copyfile('data/genes.fasta', tmp_fasta)

        # Read in original file data
        ref_seq_map = OrderedDict()
        with Fasta('data/genes.fasta', as_raw=True, strict_bounds=True) as fasta:
            for name, seq in fasta.records.items():
                ref_seq_map[name] = seq[:].lower()

        # Now write file, using many threads and simultaneously reading each sequence in pieces
        with Fasta(tmp_fasta, as_raw=True, strict_bounds=True, mutable=True) as fasta:
            threads = []
            for name, seq in fasta.records.items():
                threads.append(_ThreadWriteSequence(self.rand, name, seq))

            for thread in threads:
                thread.start()

            for thread in threads:
                thread.join()

            fasta.faidx.file.flush()

        # Now read written Fasta file, and compare it to the original
        thread_write_seq_map = OrderedDict()
        with Fasta(tmp_fasta, as_raw=True, strict_bounds=True) as fasta:
            for name, seq in fasta.records.items():
                thread_write_seq_map[name] = seq[:]

        self.assertEqual(thread_write_seq_map, ref_seq_map)

    def test_simultaneous_reads_and_writes(self):
        """
        Combine the above two tests to check that interleaved reads and writes don't conflict.
        """

        tmp_dir = self.tmp_dir

        tmp_fasta = os.path.join(tmp_dir, 'genes_write.fasta')
        shutil.copyfile('data/genes.fasta', tmp_fasta)

        # Read in original file data
        ref_seq_map = OrderedDict()
        with Fasta('data/genes.fasta', as_raw=True, strict_bounds=True) as fasta:
            for name, seq in fasta.records.items():
                ref_seq_map[name] = seq[:].lower()

        # Initialize map with fasta sequence names to enforce same ordering as 'ref_seq_map'
        thread_result_lock = threading.Lock()
        thread_read_seq_map = OrderedDict((name, None) for name in ref_seq_map)

        # Now write file, using many threads and simultaneously reading each sequence in pieces
        with Fasta(tmp_fasta, as_raw=True, strict_bounds=True, mutable=True) as fasta:
            threads = []
            for name, seq in fasta.records.items():
                threads.append(_ThreadWriteSequence(self.rand, name, seq))
                threads.append(_ThreadReadSequence(self.rand, thread_read_seq_map, thread_result_lock, name, seq))

            for thread in threads:
                thread.start()

            for thread in threads:
                thread.join()

            fasta.faidx.file.flush()

        # Now read written Fasta file, and compare it to the original
        thread_write_seq_map = OrderedDict()
        with Fasta(tmp_fasta, as_raw=True, strict_bounds=True) as fasta:
            for name, seq in fasta.records.items():
                thread_write_seq_map[name] = seq[:]

        # Change read strings to lower-case (may be a mixture of lower and upper)
        for name in ref_seq_map.keys():
            thread_read_seq_map[name] = thread_read_seq_map[name].lower()

        self.assertEqual(thread_write_seq_map, ref_seq_map)
        self.assertEqual(thread_read_seq_map, ref_seq_map)

