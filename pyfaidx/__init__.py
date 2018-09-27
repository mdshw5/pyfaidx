# pylint: disable=R0913, R0914, C0301
"""
Fasta file -> Faidx -> Fasta -> FastaRecord -> Sequence
"""

from __future__ import division
import os
import sys
from os.path import getmtime
from six import PY2, PY3, string_types, integer_types
from six.moves import zip_longest
try:
    from collections import OrderedDict
except ImportError:  #python 2.6
    from ordereddict import OrderedDict
from collections import namedtuple
import re
import string
import warnings
from math import ceil
from threading import Lock
from smart_open import smart_open

if sys.version_info > (3, ):
    buffer = memoryview

dna_bases = re.compile(r'([ACTGNactgnYRWSKMDVHBXyrwskmdvhbx]+)')

__version__ = '0.5.5'


class KeyFunctionError(ValueError):
    """Raised if the key_function argument is invalid."""


class FastaIndexingError(Exception):
    """Raised if we encounter malformed FASTA that prevents indexing."""


class IndexNotFoundError(IOError):
    """Raised if read_fai cannot open the index file."""


class FastaNotFoundError(IOError):
    """Raised if the fasta file cannot be opened."""


class FetchError(IndexError):
    """Raised if a request to fetch a FASTA sequence cannot be fulfilled."""


class BedError(ValueError):
    """Indicates a malformed BED entry."""


class RegionError(Exception):
    # This exception class is currently unused, but has been retained for
    # backwards compatibility.
    """A region error occurred."""


class UnsupportedCompressionFormat(IOError):
    """
    Raised when a FASTA file is given with a recognized but unsupported
    compression extension.
    """


class Sequence(object):
    """
    name = FASTA entry name
    seq = FASTA sequence
    start, end = coordinates of subsequence (optional)
    comp = boolean switch for complement property
    """

    def __init__(self, name='', seq='', start=None, end=None, comp=False):
        self.name = name
        self.seq = seq
        self.start = start
        self.end = end
        self.comp = comp
        assert isinstance(name, string_types)
        assert isinstance(seq, string_types)

    def __getitem__(self, n):
        """ Returns a sliced version of Sequence
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=1, end=6)
        >>> x
        >chr1:1-6
        ATCGTA
        >>> x[:3]
        >chr1:1-3
        ATC
        >>> x[3:]
        >chr1:4-6
        GTA
        >>> x[1:-1]
        >chr1:2-5
        TCGT
        >>> x[::-1]
        >chr1:6-1
        ATGCTA
        >>> x[::-3]
        >chr1
        AC
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=0, end=6)
        >>> x
        >chr1:0-6
        ATCGTA
        >>> x[:3]
        >chr1:0-3
        ATC
        >>> x[3:]
        >chr1:3-6
        GTA
        >>> x[1:-1]
        >chr1:1-5
        TCGT
        >>> x[::-1]
        >chr1:6-0
        ATGCTA
        >>> x[::-3]
        >chr1
        AC
        """
        if self.start is None or self.end is None:
            correction_factor = 0
        elif len(
                self.seq
        ) == abs(self.end - self.start) + 1:  # determine coordinate system
            one_based = True
            correction_factor = -1
        elif len(self.seq) == abs(self.end - self.start):
            one_based = False
            correction_factor = 0
        elif len(self.seq) != abs(self.end - self.start):
            raise ValueError(
                "Coordinates (Sequence.start=%s and Sequence.end=%s) imply a different length than Sequence.seq (len=%s). Did you modify Sequence.seq?"
                % (self.start, self.end, len(self.seq)))

        if isinstance(n, slice):
            slice_start, slice_stop, slice_step = n.indices(len(self))
            if self.start is None or self.end is None:  # there should never be self.start != self.end == None
                start = None
                end = None
                return self.__class__(self.name, self.seq[n], start, end,
                                      self.comp)
            self_end, self_start = (self.end, self.start)
            if abs(slice_step) > 1:
                start = None
                end = None
            elif slice_step == -1:  # flip the coordinates when we reverse
                if slice_stop == -1:
                    slice_stop = 0
                start = self_end - slice_stop
                end = self_start + slice_stop
                #print(locals())
            else:
                start = self_start + slice_start
                end = self_start + slice_stop + correction_factor
            return self.__class__(self.name, self.seq[n], start, end,
                                  self.comp)
        elif isinstance(n, integer_types):
            if n < 0:
                n = len(self) + n
            if self.start:
                return self.__class__(self.name, self.seq[n], self.start + n,
                                      self.start + n, self.comp)
            else:
                return self.__class__(self.name, self.seq[n], self.comp)

    def __str__(self):
        return self.seq

    def __neg__(self):
        """ Returns the reverse compliment of sequence
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=1, end=6)
        >>> x
        >chr1:1-6
        ATCGTA
        >>> y = -x
        >>> y
        >chr1:6-1 (complement)
        TACGAT
        >>> -y
        >chr1:1-6
        ATCGTA
        """
        return self[::-1].complement

    def __repr__(self):
        return '\n'.join([''.join(['>', self.fancy_name]), self.seq])

    def __len__(self):
        """
        >>> len(Sequence('chr1', 'ACT'))
        3
        """
        return len(self.seq)

    @property
    def fancy_name(self):
        """ Return the fancy name for the sequence, including start, end, and complementation.
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=1, end=6, comp=True)
        >>> x.fancy_name
        'chr1:1-6 (complement)'
        """
        name = self.name
        if self.start is not None and self.end is not None:
            name = ':'.join([name, '-'.join([str(self.start), str(self.end)])])
        if self.comp:
            name += ' (complement)'
        return name

    @property
    def long_name(self):
        """ DEPRECATED: Use fancy_name instead.
        Return the fancy name for the sequence, including start, end, and complementation.
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=1, end=6, comp=True)
        >>> x.long_name
        'chr1:1-6 (complement)'
        """
        msg = "The `Sequence.long_name` property is deprecated, and will be removed in future versions. Please use `Sequence.fancy_name` instead."
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        return self.fancy_name

    @property
    def complement(self):
        """ Returns the compliment of self.
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> x.complement
        >chr1 (complement)
        TAGCAT
        """
        comp = self.__class__(
            self.name, complement(self.seq), start=self.start, end=self.end)
        comp.comp = False if self.comp else True
        return comp

    @property
    def reverse(self):
        """ Returns the reverse of self.
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> x.reverse
        >chr1
        ATGCTA
        """
        return self[::-1]

    @property
    def orientation(self):
        """ get the orientation forward=1, reverse=-1
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=1, end=6)
        >>> x.orientation
        1
        >>> x.complement.orientation is None
        True
        >>> x[::-1].orientation is None
        True
        >>> x = -x
        >>> x.orientation
        -1
        """
        if self.start < self.end and not self.comp:
            return 1
        elif self.start > self.end and self.comp:
            return -1
        else:
            return None

    @property
    def gc(self):
        """ Return the GC content of seq as a float
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> y = round(x.gc, 2)
        >>> y == 0.33
        True
        """
        g = self.seq.count('G')
        g += self.seq.count('g')
        c = self.seq.count('C')
        c += self.seq.count('c')
        return (g + c) / len(self.seq)


class IndexRecord(
        namedtuple('IndexRecord',
                   ['rlen', 'offset', 'lenc', 'lenb', 'bend', 'prev_bend'])):
    __slots__ = ()

    def __getitem__(self, key):
        if type(key) == str:
            return getattr(self, key)
        return tuple.__getitem__(self, key)

    def __str__(self):
        return "{rlen:n}\t{offset:n}\t{lenc:n}\t{lenb:n}\n".format(
            **self._asdict())

    def __len__(self):
        return self.rlen


class Faidx(object):
    """ A python implementation of samtools faidx FASTA indexing """

    def __init__(self,
                 filename,
                 default_seq=None,
                 key_function=lambda x: x,
                 as_raw=False,
                 strict_bounds=False,
                 read_ahead=None,
                 mutable=False,
                 split_char=None,
                 duplicate_action="stop",
                 filt_function=lambda x: True,
                 one_based_attributes=True,
                 read_long_names=False,
                 sequence_always_upper=False,
                 rebuild=True,
                 build_index=True):
        """
        filename: name of fasta file
        key_function: optional callback function which should return a unique
          key for the self.index dictionary when given rname.
        as_raw: optional parameter to specify whether to return sequences as a
          Sequence() object or as a raw string.
          Default: False (i.e. return a Sequence() object).
        """
        self.filename = filename

        if filename.lower().endswith('.bgz') or filename.lower().endswith(
                '.gz'):
            # Only try to import Bio if we actually need the bgzf reader.
            try:
                from Bio import bgzf
                from Bio import __version__ as bgzf_version
                from distutils.version import LooseVersion
                if LooseVersion(bgzf_version) < LooseVersion('1.73'):
                    raise ImportError
            except ImportError:
                raise ImportError(
                    "BioPython >= 1.73 must be installed to read block gzip files.")
            else:
                self._fasta_opener = bgzf.open
                self._bgzf = True
        elif filename.lower().endswith('.bz2') or filename.lower().endswith(
                '.zip'):
            raise UnsupportedCompressionFormat(
                "Compressed FASTA is only supported in BGZF format. Use "
                "bgzip to compresss your FASTA.")
        else:
            self._fasta_opener = smart_open
            self._bgzf = False

        try:
            self.file = self._fasta_opener(filename, 'r+b'
                                           if mutable else 'rb')
        except (ValueError, IOError) as e:
            if str(e).find('BGZF') > -1:
                raise UnsupportedCompressionFormat(
                    "Compressed FASTA is only supported in BGZF format. Use "
                    "the samtools bgzip utility (instead of gzip) to "
                    "compress your FASTA.")
            else:
                raise FastaNotFoundError(
                    "Cannot read FASTA file %s" % filename)

        self.indexname = filename + '.fai'
        self.read_long_names = read_long_names
        self.key_function = key_function
        try:
            key_fn_test = self.key_function(
                "TestingReturnType of_key_function")
            if not isinstance(key_fn_test, string_types):
                raise KeyFunctionError(
                    "key_function argument should return a string, not {0}".
                    format(type(key_fn_test)))
        except Exception as e:
            pass
        self.filt_function = filt_function
        assert duplicate_action in ("stop", "first", "last", "longest",
                                    "shortest", "drop")
        self.duplicate_action = duplicate_action
        self.as_raw = as_raw
        self.default_seq = default_seq
        if self._bgzf and self.default_seq is not None:
            raise FetchError(
                "The default_seq argument is not supported with using BGZF compression. Please decompress your FASTA file and try again."
            )
        if self._bgzf:
            self.strict_bounds = True
        else:
            self.strict_bounds = strict_bounds
        self.split_char = split_char
        self.one_based_attributes = one_based_attributes
        self.sequence_always_upper = sequence_always_upper
        self.index = OrderedDict()
        self.lock = Lock()
        self.buffer = dict((('seq', None), ('name', None), ('start', None),
                            ('end', None)))
        if not read_ahead or isinstance(read_ahead, integer_types):
            self.read_ahead = read_ahead
        elif not isinstance(read_ahead, integer_types):
            raise ValueError("read_ahead value must be int, not {0}".format(
                type(read_ahead)))

        self.mutable = mutable
        with self.lock:  # lock around index generation so only one thread calls method
            try:
                if os.path.exists(self.indexname) and getmtime(
                        self.indexname) >= getmtime(self.filename):
                    self.read_fai()
                elif os.path.exists(self.indexname) and getmtime(
                        self.indexname) < getmtime(
                            self.filename) and not rebuild:
                    self.read_fai()
                    warnings.warn(
                        "Index file {0} is older than FASTA file {1}.".format(
                            self.indexname, self.filename), RuntimeWarning)
                elif build_index:
                    self.build_index()
                    self.read_fai()
                else:
                    self.read_fai()

            except FastaIndexingError:
                os.remove(self.indexname)
                self.file.close()
                raise
            except Exception:
                # Handle potential exceptions other than 'FastaIndexingError'
                self.file.close()
                raise

    def __contains__(self, region):
        if not self.buffer['name']:
            return False
        name, start, end = region
        if self.buffer['name'] == name and self.buffer['start'] <= start and self.buffer['end'] >= end:
            return True
        else:
            return False

    def __repr__(self):
        return 'Faidx("%s")' % (self.filename)

    def _index_as_string(self):
        """ Returns the string representation of the index as iterable """
        for k, v in self.index.items():
            yield '\t'.join([k, str(v)])

    def read_fai(self):
        try:
            with smart_open(self.indexname, encoding='utf8') as index:
                prev_bend = 0
                drop_keys = []
                for line in index:
                    line = line.rstrip()
                    rname, rlen, offset, lenc, lenb = line.split('\t')
                    rlen, offset, lenc, lenb = map(int,
                                                   (rlen, offset, lenc, lenb))
                    newlines = int(ceil(rlen / lenc) * (lenb - lenc))
                    bend = offset + newlines + rlen
                    rec = IndexRecord(rlen, offset, lenc, lenb, bend,
                                      prev_bend)
                    if self.read_long_names:
                        rname = self._long_name_from_index_record(rec)
                    if self.split_char:
                        rname = filter(self.filt_function,
                                       self.key_function(rname).split(
                                           self.split_char))
                    else:
                        # filter must act on an iterable
                        rname = filter(self.filt_function,
                                       [self.key_function(rname)])
                    for key in rname:  # mdshw5/pyfaidx/issues/64
                        if key in self.index:
                            if self.duplicate_action == "stop":
                                raise ValueError('Duplicate key "%s"' % key)
                            elif self.duplicate_action == "first":
                                continue
                            elif self.duplicate_action == "last":
                                self.index[key] = rec
                            elif self.duplicate_action == "longest":
                                if len(rec) > len(self.index[key]):
                                    self.index[key] = rec
                            elif self.duplicate_action == "shortest":
                                if len(rec) < len(self.index[key]):
                                    self.index[key] = rec
                            elif self.duplicate_action == "drop":
                                if key not in drop_keys:
                                    drop_keys.append(key)
                        else:
                            self.index[key] = rec
                    prev_bend = bend
            for dup in drop_keys:
                self.index.pop(dup, None)
        except IOError:
            raise IndexNotFoundError(
                "Could not read index file %s" % self.indexname)

    def build_index(self):
        try:
            with self._fasta_opener(self.filename, 'rb') as fastafile:
                with smart_open(self.indexname, 'w') as indexfile:
                    rname = None  # reference sequence name
                    offset = 0  # binary offset of end of current line
                    rlen = 0  # reference character length
                    blen = None  # binary line length (includes newline)
                    clen = None  # character line length
                    bad_lines = []  # lines > || < than blen
                    thisoffset = offset

                    lastline = None
                    for i, line in enumerate(fastafile):
                        line_blen = len(line)
                        line = line.decode()
                        line_clen = len(line.rstrip('\n\r'))
                        lastline = i
                        # write an index line
                        if line[0] == '>':
                            valid_entry = check_bad_lines(
                                rname, bad_lines, i - 1)
                            if valid_entry and i > 0:
                                indexfile.write(
                                    "{0}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\n".format(
                                        rname, rlen, thisoffset, clen, blen))
                            elif not valid_entry:
                                raise FastaIndexingError(
                                    "Line length of fasta"
                                    " file is not "
                                    "consistent! "
                                    "Inconsistent line found in >{0} at "
                                    "line {1:n}.".format(
                                        rname, bad_lines[0][0] + 1))
                            blen = None
                            rlen = 0
                            clen = None
                            bad_lines = []
                            try:  # must catch empty deflines (actually these might be okay: https://github.com/samtools/htslib/pull/258)
                                rname = line.rstrip('\n\r')[1:].split()[
                                    0]  # duplicates are detected with read_fai
                            except IndexError:
                                raise FastaIndexingError(
                                    "Bad sequence name %s at line %s." %
                                    (line.rstrip('\n\r'), str(i)))
                            offset += line_blen
                            thisoffset = fastafile.tell(
                            ) if self._bgzf else offset
                        else:  # check line and advance offset
                            if not blen:
                                blen = line_blen
                            if not clen:
                                clen = line_clen
                            # only one short line should be allowed
                            # before we hit the next header, and it
                            # should be the last line in the entry
                            if line_blen != blen or line_blen == 1:
                                bad_lines.append((i, line_blen))
                            offset += line_blen
                            rlen += line_clen

                    # write the final index line, if there is one.
                    if lastline is not None:
                        valid_entry = check_bad_lines(
                            rname, bad_lines, lastline
                        )  # advance index since we're at the end of the file
                        if valid_entry:
                            indexfile.write(
                                "{0:s}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\n".format(
                                    rname, rlen, thisoffset, clen, blen))
                        else:
                            raise FastaIndexingError(
                                "Line length of fasta"
                                " file is not "
                                "consistent! "
                                "Inconsistent line found in >{0} at "
                                "line {1:n}.".format(rname,
                                                     bad_lines[0][0] + 1))
        except (IOError, FastaIndexingError) as e:
            if isinstance(e, IOError):
                raise IOError(
                    "%s may not be writable. Please use Fasta(rebuild=False), Faidx(rebuild=False) or faidx --no-rebuild."
                    % self.indexname)
            elif isinstance(e, FastaIndexingError):
                raise e

    def write_fai(self):
        with self.lock:
            with smart_open(self.indexname, 'w') as outfile:
                for line in self._index_as_string:
                    outfile.write(line)

    def from_buffer(self, start, end):
        i_start = start - self.buffer['start']  # want [0, 1) coordinates from [1, 1] coordinates
        i_end = end - self.buffer['start'] + 1
        return self.buffer['seq'][i_start:i_end]

    def fill_buffer(self, name, start, end):
        try:
            seq = self.from_file(name, start, end)
            self.buffer['seq'] = seq
            self.buffer['start'] = start
            self.buffer['end'] = end
            self.buffer['name'] = name
        except FetchError:
            pass

    def fetch(self, name, start, end):
        if self.read_ahead and not (name, start, end) in self:
            self.fill_buffer(name, start, end + self.read_ahead)

        if (name, start, end) in self:
            seq = self.from_buffer(start, end)
        else:
            seq = self.from_file(name, start, end)

        return self.format_seq(seq, name, start, end)

    def from_file(self, rname, start, end, internals=False):
        """ Fetch the sequence ``[start:end]`` from ``rname`` using 1-based coordinates
        1. Count newlines before start
        2. Count newlines to end
        3. Difference of 1 and 2 is number of newlines in [start:end]
        4. Seek to start position, taking newlines into account
        5. Read to end position, return sequence
        """
        assert start == int(start)
        assert end == int(end)
        try:
            i = self.index[rname]
        except KeyError:
            raise FetchError("Requested rname {0} does not exist! "
                             "Please check your FASTA file.".format(rname))
        start0 = start - 1  # make coordinates [0,1)
        if start0 < 0:
            raise FetchError(
                "Requested start coordinate must be greater than 1.")
        seq_len = end - start0

        # Calculate offset (https://github.com/samtools/htslib/blob/20238f354894775ed22156cdd077bc0d544fa933/faidx.c#L398)
        newlines_before = int(
            (start0 - 1) / i.lenc * (i.lenb - i.lenc)) if start0 > 0 else 0
        newlines_to_end = int(end / i.lenc * (i.lenb - i.lenc))
        newlines_inside = newlines_to_end - newlines_before
        seq_blen = newlines_inside + seq_len
        bstart = i.offset + newlines_before + start0
        if seq_blen < 0 and self.strict_bounds:
            raise FetchError("Requested coordinates start={0:n} end={1:n} are "
                             "invalid.\n".format(start, end))
        elif end > i.rlen and self.strict_bounds:
            raise FetchError("Requested end coordinate {0:n} outside of {1}. "
                             "\n".format(end, rname))

        with self.lock:
            if self._bgzf:  # We can't add to virtual offsets, so we need to read from the beginning of the record and trim the beginning if needed
                self.file.seek(i.offset)
                chunk = start0 + newlines_before + newlines_inside + seq_len
                chunk_seq = self.file.read(chunk).decode()
                seq = chunk_seq[start0 + newlines_before:]
            else:
                self.file.seek(bstart)

                if bstart + seq_blen > i.bend and not self.strict_bounds:
                    seq_blen = i.bend - bstart

                if seq_blen > 0:
                    seq = self.file.read(seq_blen).decode()
                elif seq_blen <= 0 and not self.strict_bounds:
                    seq = ''

        if not internals:
            return seq.replace('\n', '').replace('\r', '')
        else:
            return (seq, locals())

    def format_seq(self, seq, rname, start, end):
        start0 = start - 1
        if len(
                seq
        ) < end - start0 and self.default_seq:  # Pad missing positions with default_seq
            pad_len = end - start0 - len(seq)
            seq = ''.join([seq, pad_len * self.default_seq])
        else:  # Return less than requested range
            end = start0 + len(seq)

        if self.sequence_always_upper:
            seq = seq.upper()

        if not self.one_based_attributes:
            start = start0

        if self.as_raw:
            return seq
        else:
            return Sequence(
                name=rname, start=int(start), end=int(end), seq=seq)

    def to_file(self, rname, start, end, seq):
        """ Write sequence in region from start-end, overwriting current
        contents of the FASTA file. """
        if not self.mutable:
            raise IOError(
                "Write attempted for immutable Faidx instance. Set mutable=True to modify original FASTA."
            )
        file_seq, internals = self.from_file(rname, start, end, internals=True)

        with self.lock:
            if len(seq) != len(file_seq) - internals['newlines_inside']:
                raise IOError(
                    "Specified replacement sequence needs to have the same length as original."
                )
            elif len(seq) == len(file_seq) - internals['newlines_inside']:
                line_len = internals['i'].lenc
                if '\r\n' in file_seq:
                    newline_char = '\r\n'
                elif '\r' in file_seq:
                    newline_char = '\r'
                else:
                    newline_char = '\n'
                self.file.seek(internals['bstart'])
                if internals['newlines_inside'] == 0:
                    self.file.write(seq.encode())
                elif internals['newlines_inside'] > 0:
                    n = 0
                    m = file_seq.index(newline_char)
                    while m < len(seq):
                        self.file.write(''.join([seq[n:m], newline_char]).encode())
                        n = m
                        m += line_len
                    self.file.write(seq[n:].encode())
                    self.file.flush()

    def get_long_name(self, rname):
        """ Return the full sequence defline and description. External method using the self.index """
        index_record = self.index[rname]
        if self._bgzf:
            return self._long_name_from_bgzf(index_record)
        else:
            return self._long_name_from_index_record(index_record)

    def _long_name_from_index_record(self, index_record):
        """ Return the full sequence defline and description. Internal method passing IndexRecord """
        prev_bend = index_record.prev_bend
        defline_end = index_record.offset
        self.file.seek(prev_bend)
        return self.file.read(defline_end - prev_bend).decode()[1:-1]

    def _long_name_from_bgzf(self, index_record):
        """ Return the full sequence defline and description. Internal method passing IndexRecord
        This method is present for compatibility with BGZF files, since we cannot subtract their offsets.
        It may be possible to implement a more efficient method. """
        raise NotImplementedError(
            "FastaRecord.long_name and Fasta(read_long_names=True) "
            "are not supported currently for BGZF compressed files.")
        prev_bend = index_record.prev_bend
        self.file.seek(prev_bend)
        defline = []
        while True:
            chunk = self.file.read(4096).decode()
            defline.append(chunk)
            if '\n' in chunk or '\r' in chunk:
                break
        return ''.join(defline)[1:].split('\n\r')[0]

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class FastaRecord(object):
    __slots__ = ['name', '_fa']

    def __init__(self, name, fa):
        self.name = name
        self._fa = fa

    def __getitem__(self, n):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        try:
            if isinstance(n, slice):
                start, stop, step = n.start, n.stop, n.step
                if start is None:
                    start = 0
                if stop is None:
                    stop = len(self)
                if stop < 0:
                    stop = len(self) + stop
                if start < 0:
                    start = len(self) + start
                return self._fa.get_seq(self.name, start + 1, stop)[::step]

            elif isinstance(n, integer_types):
                if n < 0:
                    n = len(self) + n
                return self._fa.get_seq(self.name, n + 1, n + 1)
        except FetchError:
            raise

    def __iter__(self):
        """ Construct a line-based iterator that respects the original line lengths. """
        line_len = self._fa.faidx.index[self.name].lenc
        start = 0
        while True:
            end = start + line_len
            if end < len(self):
                yield self[start:end]
            else:
                yield self[start:]
                raise StopIteration
            start += line_len

    def __reversed__(self):
        """ Reverse line-based iterator """
        line_len = self._fa.faidx.index[self.name].lenc
        # have to determine last line length
        last_line = len(self) % line_len
        if last_line == 0:
            last_line = line_len
        end = len(self)
        start = end - last_line
        while True:
            if start > 0:
                yield self[start:end][::-1]
            else:
                yield self[:end][::-1]
                raise StopIteration
            if end == len(self):  # first iteration
                end -= last_line
            else:
                end -= line_len
            start = end - line_len

    def __repr__(self):
        return 'FastaRecord("%s")' % (self.name)

    def __len__(self):
        return self._fa.faidx.index[self.name].rlen

    @property
    def unpadded_len(self):
        """ Returns the length of the contig without 5' and 3' N padding.
        Functions the same as contigNonNSize in Fasta.cpp at
        https://github.com/Illumina/hap.py/blob/master/src/c%2B%2B/lib/tools/Fasta.cpp#L284
        """
        length = len(self)
        stop = False
        for line in iter(self):
            if stop:
                break
            if isinstance(line, Sequence):
                line = line.seq
            for base in line.upper():
                if base == 'N':
                    length -= 1
                else:
                    stop = True
                    break
        stop = False
        for line in reversed(self):
            if stop:
                break
            if isinstance(line, Sequence):
                line = line.seq
            for base in line.upper():
                if base == 'N':
                    length -= 1
                else:
                    stop = True
                    break
        return length

    def __str__(self):
        return str(self[:])

    @property
    def variant_sites(self):
        if isinstance(self._fa, FastaVariant):
            pos = []
            var = self._fa.vcf.fetch(self.name, 0, len(self))
            for site in var:
                if site.is_snp:
                    sample = site.genotype(self._fa.sample)
                    if sample.gt_type in self._fa.gt_type and eval(
                            self._fa.filter):
                        pos.append(site.POS)
            return tuple(pos)
        else:
            raise NotImplementedError(
                "variant_sites() only valid for FastaVariant.")

    @property
    def long_name(self):
        """ Read the actual defline from self._fa.faidx mdshw5/pyfaidx#54 """
        return self._fa.faidx.get_long_name(self.name)

    @property
    def __array_interface__(self):
        """ Implement numpy array interface for issue #139"""
        return {
            'shape': (len(self), ),
            'typestr': '|S1',
            'version': 3,
            'data': buffer(str(self).encode('ascii'))
        }


class MutableFastaRecord(FastaRecord):
    def __init__(self, name, fa):
        super(MutableFastaRecord, self).__init__(name, fa)
        if self._fa.faidx._fasta_opener != smart_open:
            raise UnsupportedCompressionFormat(
                "BGZF compressed FASTA is not supported for MutableFastaRecord. "
                "Please decompress your FASTA file.")

    def __setitem__(self, n, value):
        """Mutate sequence in region [start, end)
        to value.
        Coordinates are 0-based, end-exclusive."""
        try:
            if isinstance(n, slice):
                start, stop, step = n.start, n.stop, n.step
                if step:
                    raise IndexError("Step operator is not implemented.")
                if not start:
                    start = 0
                if not stop:
                    stop = len(self)
                if stop < 0:
                    stop = len(self) + stop
                if start < 0:
                    start = len(self) + start
                self._fa.faidx.to_file(self.name, start + 1, stop, value)

            elif isinstance(n, integer_types):
                if n < 0:
                    n = len(self) + n
                return self._fa.faidx.to_file(self.name, n + 1, n + 1, value)
        except (FetchError, IOError):
            raise


class Fasta(object):
    def __init__(self,
                 filename,
                 default_seq=None,
                 key_function=lambda x: x,
                 as_raw=False,
                 strict_bounds=False,
                 read_ahead=None,
                 mutable=False,
                 split_char=None,
                 filt_function=lambda x: True,
                 one_based_attributes=True,
                 read_long_names=False,
                 duplicate_action="stop",
                 sequence_always_upper=False,
                 rebuild=True,
                 build_index=True):
        """
        An object that provides a pygr compatible interface.
        filename: name of fasta file
        """
        self.filename = filename
        self.mutable = mutable
        self.faidx = Faidx(
            filename,
            key_function=key_function,
            as_raw=as_raw,
            default_seq=default_seq,
            strict_bounds=strict_bounds,
            read_ahead=read_ahead,
            mutable=mutable,
            split_char=split_char,
            filt_function=filt_function,
            one_based_attributes=one_based_attributes,
            read_long_names=read_long_names,
            duplicate_action=duplicate_action,
            sequence_always_upper=sequence_always_upper,
            rebuild=rebuild,
            build_index=build_index)
        self.keys = self.faidx.index.keys
        if not self.mutable:
            self.records = dict(
                [(rname, FastaRecord(rname, self)) for rname in self.keys()])
        elif self.mutable:
            self.records = dict([(rname, MutableFastaRecord(rname, self))
                                 for rname in self.keys()])

    def __contains__(self, rname):
        """Return True if genome contains record."""
        return rname in self.faidx.index

    def __getitem__(self, rname):
        """Return a chromosome by its name, or its numerical index."""
        if isinstance(rname, integer_types):
            rname = tuple(self.keys())[rname]
        try:
            return self.records[rname]
        except KeyError:
            raise KeyError("{0} not in {1}.".format(rname, self.filename))

    def __repr__(self):
        return 'Fasta("%s")' % (self.filename)

    def __iter__(self):
        for rname in self.keys():
            yield self[rname]

    def get_seq(self, name, start, end, rc=False):
        """Return a sequence by record name and interval [start, end).

        Coordinates are 1-based, end-exclusive.
        If rc is set, reverse complement will be returned.
        """
        # Get sequence from real genome object and save result.
        seq = self.faidx.fetch(name, start, end)
        if rc:
            return -seq
        else:
            return seq

    def get_spliced_seq(self, name, intervals, rc=False):
        """Return a sequence by record name and list of intervals

        Interval list is an iterable of [start, end].
        Coordinates are 1-based, end-exclusive.
        If rc is set, reverse complement will be returned.
        """
        # Get sequence for all intervals
        chunks = [self.faidx.fetch(name, s, e) for s, e in intervals]
        start = chunks[0].start
        end = chunks[-1].end

        # reverce complement
        if rc:
            seq = "".join([(-chunk).seq for chunk in chunks[::-1]])
        else:
            seq = "".join([chunk.seq for chunk in chunks])

        # Sequence coordinate validation wont work since
        # len(Sequence.seq) != end - start
        return Sequence(name=name, seq=seq, start=None, end=None)

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.faidx.__exit__(*args)


class FastaVariant(Fasta):
    """ Return consensus sequence from FASTA and VCF inputs
    """
    expr = set(('>', '<', '=', '!'))

    def __init__(self,
                 filename,
                 vcf_file,
                 sample=None,
                 het=True,
                 hom=True,
                 call_filter=None,
                 **kwargs):
        super(FastaVariant, self).__init__(filename, **kwargs)
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam must be installed for FastaVariant.")
        try:
            import vcf
        except ImportError:
            raise ImportError("PyVCF must be installed for FastaVariant.")
        if call_filter is not None:
            try:
                key, expr, value = call_filter.split()  # 'GQ > 30'
            except IndexError:
                raise ValueError(
                    "call_filter must be a string in the format 'XX <>!= NN'")
            assert all([x in self.expr for x in list(expr)])
            assert all([x in string.ascii_uppercase for x in list(key)])
            assert all([x in string.printable for x in list(value)])
            self.filter = "sample['{key}'] {expr} {value}".format(**locals())
        else:
            self.filter = 'True'
        if os.path.exists(vcf_file):
            self.vcf = vcf.Reader(filename=vcf_file)
        else:
            raise IOError("File {0} does not exist.".format(vcf_file))
        if sample is not None:
            self.sample = sample
        else:
            self.sample = self.vcf.samples[0]
            if len(self.vcf.samples) > 1 and sample is None:
                warnings.warn("Using sample {0} genotypes.".format(
                    self.sample), RuntimeWarning)
        if het and hom:
            self.gt_type = set((1, 2))
        elif het:
            self.gt_type = set((1, ))
        elif hom:
            self.gt_type = set((2, ))
        else:
            self.gt_type = set()

    def __repr__(self):
        return 'FastaVariant("%s", "%s", gt="%s")' % (self.filename,
                                                      self.vcf.filename,
                                                      str(self.gt_type))

    def get_seq(self, name, start, end):
        """Return a sequence by record name and interval [start, end).
        Replace positions with polymorphism with variant.
        Coordinates are 0-based, end-exclusive.
        """
        seq = self.faidx.fetch(name, start, end)
        if self.faidx.as_raw:
            seq_mut = list(seq)
            del seq
        else:
            seq_mut = list(seq.seq)
            del seq.seq
        var = self.vcf.fetch(name, start - 1, end)
        for record in var:
            if record.is_snp:  # skip indels
                sample = record.genotype(self.sample)
                if sample.gt_type in self.gt_type and eval(self.filter):
                    alt = record.ALT[0]
                    i = (record.POS - 1) - (start - 1)
                    seq_mut[i:i + len(alt)] = str(alt)
        # slice the list in case we added an MNP in last position
        if self.faidx.as_raw:
            return ''.join(seq_mut[:end - start + 1])
        else:
            seq.seq = ''.join(seq_mut[:end - start + 1])
            return seq


def wrap_sequence(n, sequence, fillvalue=''):
    args = [iter(sequence)] * n
    for line in zip_longest(fillvalue=fillvalue, *args):
        yield ''.join(line + ("\n", ))


# To take a complement, we map each character in the first string in this pair
# to the corresponding character in the second string.
complement_map = ('ACTGNactgnYRWSKMDVHBXyrwskmdvhbx',
                  'TGACNtgacnRYWSMKHBDVXrywsmkhbdvx')
invalid_characters_set = set(
    chr(x) for x in range(256) if chr(x) not in complement_map[0])
invalid_characters_string = ''.join(invalid_characters_set)

if PY3:
    complement_table = str.maketrans(complement_map[0], complement_map[1],
                                     invalid_characters_string)
    translate_arguments = (complement_table, )
elif PY2:
    complement_table = string.maketrans(complement_map[0], complement_map[1])
    translate_arguments = (complement_table, invalid_characters_string)


def complement(seq):
    """ Returns the complement of seq.
    >>> seq = 'ATCGTA'
    >>> complement(seq)
    'TAGCAT'
    """
    seq = str(seq)
    result = seq.translate(*translate_arguments)
    if len(result) != len(seq):
        first_invalid_position = next(
            i for i in range(len(seq)) if seq[i] in invalid_characters_set)
        raise ValueError(
            "Sequence contains non-DNA character '{0}' at position {1:n}\n".
            format(seq[first_invalid_position], first_invalid_position + 1))
    return result


def translate_chr_name(from_name, to_name):
    chr_name_map = dict(zip(from_name, to_name))

    def map_to_function(rname):
        return chr_name_map[rname]

    return map_to_function


def bed_split(bed_entry):
    try:
        rname, start, end = bed_entry.rstrip().split()[:3]
    except (IndexError, ValueError):
        raise BedError('Malformed BED entry! {0}\n'.format(bed_entry.rstrip()))
    start, end = (int(start), int(end))
    return (rname, start, end)


def ucsc_split(region):
    try:
        rname, interval = region.split(':')
    except ValueError:
        rname = region
        interval = None
    try:
        start, end = interval.split('-')
        start, end = (int(start) - 1, int(end))
    except (AttributeError, ValueError):
        start, end = (None, None)
    return (rname, start, end)


def check_bad_lines(rname, bad_lines, i):
    """ Find inconsistent line lengths in the middle of an
    entry. Allow blank lines between entries, and short lines
    occurring at the last line of an entry. Returns boolean
    validating the entry.
    >>> check_bad_lines('chr0', [(10, 79)], 10)
    True
    >>> check_bad_lines('chr0', [(9, 79)], 10)
    False
    >>> check_bad_lines('chr0', [(9, 79), (10, 1)], 10)
    True
    """
    if len(bad_lines) == 0:
        return True
    elif len(bad_lines) == 1:
        if bad_lines[0][0] == i:  # must be last line
            return True
        else:
            return False
    elif len(bad_lines) == 2:
        if bad_lines[0][0] == i:  # must not be last line
            return False
        elif bad_lines[1][0] == i and bad_lines[1][1] == 1:  # blank last line
            if bad_lines[0][0] + 1 == i and bad_lines[0][1] > 1:  # non-blank line
                return True
        else:
            return False
    if len(bad_lines) > 2:
        return False
    raise RuntimeError("Unhandled exception during fasta indexing at entry " + rname + \
                       "Please report this issue at https://github.com/mdshw5/pyfaidx/issues " + \
                       str(bad_lines))


if __name__ == "__main__":
    import doctest
    doctest.testmod()
