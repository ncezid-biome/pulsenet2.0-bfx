from collections import defaultdict, namedtuple
import os
import re
import csv
import json
from csv import DictReader
from os import getenv
from pathlib import Path
from shutil import move
from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from subprocess import run, PIPE
from itertools import zip_longest, repeat

from src.inputs import EcStxTyperInputs
from src.outputs import EcStxTyperOutputs


CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")
_QC_SECTION = "STEC"


# some tooling
_FASTAEXTS = ['.fna', '.fasta', '.fsa']

def fasta_iterator(flname):
    with open(flname, 'r') as f:
        sequence_parts = []
        key = ''
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                if key:
                    full_seqence = ''.join(sequence_parts).upper()
                    yield (key, full_seqence)
                # Start of the file
                key = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if key:
            # The last sequence in the file
            full_seqence = ''.join(sequence_parts).upper()
            yield (key, full_seqence)

def is_fasta(path):
    # Check if the file is a fasta format
    if not os.path.exists(path):
        return False
    file_extension = os.path.splitext(path)[1]
    return file_extension in _FASTAEXTS

def parse_fasta(flname, rename=False):
    if is_fasta(flname):
        fasta_sequences = {}
        if rename:
            for i, (name, sequence) in enumerate(fasta_iterator(flname), 1):
                new_name = 'contig_' + str(i)
                fasta_sequences[new_name] = sequence
            return fasta_sequences
        else:
            fasta_sequences.update(fasta_iterator(flname))
            return fasta_sequences
    else:
        raise RuntimeError('Requested file path is not of type fasta')

SequenceInfo = namedtuple('SequenceInfo', ['locus', 'allele', 'accession', 'sequence', 'other'])
LocusInfo = namedtuple('LocusInfo', ['locus', 'note', 'antibiotic', 'other'])

def sequence_parser(header, sequence, sep = ':'):
    # Default parser for the database files
    parts = header.split(sep)
    parts = list(map(str.strip, parts))
    while len(parts) < 4:
        parts.append('')
    return SequenceInfo(
        locus = parts[0],
        allele = parts[1],
        accession = parts[2],
        sequence = sequence,
        other = parts[3]
    )

def notes_parser(line, sep=':'):
    # Default parser for the notes files
    parts = list(map(str.strip, line.split(sep)))
    while len(parts) < 3:
        parts.append('')
    antibiotic = parts[1].replace('resistance', '')
    antibiotic = map(str.strip, antibiotic.split(','))

    return LocusInfo(
        locus = parts[0],
        note = parts[1],
        antibiotic = antibiotic,
        other = parts[2]
    )


class DbInfo(object):
    # Class that will hold the db information
    def __init__(self, dirpath, seq_parser=sequence_parser,
                 note_parser=notes_parser):

        self._notes = {}
        self._sequences = {}
        self._dirpath = dirpath
        self._separator = None
        # if dirpath is None or not check_dir(dirpath):
        #     raise RuntimeError('Invalid path provided for '
        #                        ' database fastas: {}'.format(str(dirpath)))
        self.load_database(dirpath, seq_parser, note_parser)

    def load_database(self, dirpath, seq_parser, note_parser):
        sequence_counts = defaultdict(dict)
        allele_id_template = '{}_{}'
        allele_id_template_i = '{}_{}-{}'
        self._separator = '-'
        for seq_file in os.listdir(dirpath):
            file_path = os.path.join(dirpath, seq_file)
            if not is_fasta(file_path):
                continue
            sequences = parse_fasta(file_path)
            for seq_id, sequence in sequences.items():
                seq_info = seq_parser(seq_id, sequence)
                allele_id = allele_id_template.format(seq_info.locus, seq_info.allele)
                if allele_id in sequence_counts:
                    # There are a lot of alleles with the same name
                    # new is created to make each sequence unique
                    new_id = allele_id_template_i.format(
                        seq_info.locus,
                        seq_info.allele,
                        len(sequence_counts[allele_id])
                    )
                    sequence_counts[allele_id][new_id] = True
                    self._sequences[new_id] = seq_info
                else:
                    sequence_counts[allele_id][allele_id] = True
                    self._sequences[allele_id] = seq_info

        # Load notes if they exist
        file_path = os.path.join(dirpath, 'notes.txt')
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line[0] == '#':
                        continue
                    notes_info = note_parser(line)
                    self._notes[notes_info.locus] = notes_info

    @property
    def sequences(self):
        return self._sequences

    @property
    def notes(self):
        return self._notes

    def export_sequences(self, filepath):
        # valid_dir(os.path.dirname(filepath))
        with open(filepath, 'w') as f:
            for seq_id, seq_info in self._sequences.items():
                ostr = '>{}|{}\n{}\n'
                f.write(ostr.format(
                    seq_id,
                    len(seq_info.sequence),
                    seq_info.sequence
                ))

    def get_refseq(self, ref):
        if ref in self._sequences:
            return self._sequences[ref].sequence
        raise KeyError('Missing reference: {}'.format(ref))

    def results_parser(self, results, f=None):
        if f is not None and callable(f):
            return f(self, results)
        sequences = self.sequences
        notes = self.notes
        if not len(sequences):
            raise RuntimeError('Trying to assimilate results for'
                               ' something with no sequences')
        results_out = {
            'results': {},
            'extra': []
        }
        for result, geno_regions in results.items():
            sequence_info = sequences[result]
            locus = sequence_info.locus
            allele = sequence_info.allele
            gene_name = locus
            results_out['results'][gene_name] = True

            hit_information = [
                {
                    'locus': locus,
                    'identity': geno_region.identity,
                    'coverage': geno_region.coverage,
                    'allele': allele,
                    'hits': [
                        {
                            'contig_id': hit.query_id,
                            'query_start': hit.query_start,
                            'query_stop': hit.query_stop,
                            'reference_start': hit.reference_start,
                            'reference_stop': hit.reference_stop,
                            'full_match': hit.full_match,
                        } for hit in geno_region.locations
                    ]
                } for geno_region in geno_regions
            ]
            results_out['extra'].extend(hit_information)
        for sequence_info in sequences.values():
            locus = sequence_info.locus
            allele = sequence_info.allele
            gene_name = locus
            if gene_name not in results_out['results']:
                results_out['results'][gene_name] = False
        return results_out


_deletions = re.compile(r'-([0-9]+)([ACGTNacgtn]+)')
_insertions = re.compile(r'\+([0-9]+)([ACGTNacgtn]+)')
_substitutions = re.compile(r'[ACGTNacgtn]')
_remove = re.compile(r'[$<>]')
_start_read = re.compile(r'(\^)(.)')
_reference = re.compile(r'[.,]')
_asterisk = re.compile(r'[\*]')

def counter(l):
    c = {}
    for i in l:
        if i in c:
            c[i] += 1
        else:
            c[i] = 1
    return c

_SET_TO_NON_ACTG = {}
def get_non_iupac(f_set):
    return _SET_TO_NON_ACTG.get(f_set, None)


class ConsensusSequence(object):
    def __init__(self):
        self._start = -1
        self._stop = -1
        self._seq = []
        self._count = -1
        self._ambiguous = {}
        self._ref = ''

    def initialize(self, ref, start, stop, nuc, count, ref_call):
        self._start = start
        self._stop = stop
        if not isinstance(nuc, ConsensusPosition):
            cp = ConsensusPosition(stop, nuc, count, ref_call)
        else:
            cp = nuc
        if cp.ambiguous:
            self._ambiguous[stop] = cp
        self._seq.append(cp)
        self._count = count
        self._ref = ref

    def add_nuc(self, ref, pos, nuc, count, ref_call):
        if self._start == -1:
            return self.initialize(ref, pos, pos, nuc, count, ref_call)
        assert self.ref == ref
        if not isinstance(nuc, ConsensusPosition):
            cp = ConsensusPosition(pos, nuc, count, ref_call)
        else:
            cp = nuc
        if not pos == self._stop + 1:
            raise RuntimeError('Adding to consensus sequence needs to be contiguous')
        if cp.ambiguous:
            self._ambiguous[pos - 1] = cp
        self._stop += 1
        self._seq.append(cp)
        self._count += count

    @staticmethod
    def merge(first, second):
        if type(first) != type(second):
            raise RuntimeError('Type mismatch: {} vs. {}'.format(
                type(first), type(second)))
        if first.start > second.start:
            return ConsensusSequence.merge(second, first)
        if first.start <= second.end <= first.end:
            raise RuntimeError('Cannot merge overlapping sequences')
        # Get the gap between the two sequences in the case
        # that we have only parts of the gene covered
        gap_count = max(0, second.start - first.end - 1)
        while gap_count:
            first.add_nuc(first.end + 1, '-', 0)
            gap_count -= 1
        # Append the sequence to the other
        for i, nuc in enumerate(second.seq, second.start):
            first.add_nuc(i, nuc, 0)
        # Add the coverage value so we can
        # get the coverage later
        first.count += second.count
        return first

    def flatten(self):
        rm = []
        for k, v in self.ambiguous.items():
            if not v.flatten():
                rm.append(k)
        for k in rm:
            del self.ambiguous[k]
        self.count = sum(p.count for p in self)

    def get_fragment(self, start=None, stop=None):
        if start is None:
            start = self.start
        if stop is None:
            stop = self.stop
        for i in range(start - 1, stop):
            if self.seq[i].ambiguous:
                if isinstance(self.seq[i].nuc, str):
                    yield self.seq[i].nuc
                else:
                    nucs_here = set(map(str.upper, self.seq[i].nuc))
                    if '-' in nucs_here:
                        yield 'N'
                    else:
                        nucs_here = set(nucs_here)
                        yield_fstr = '[{}]'
                        ins = []
                        for n in nucs_here:
                            if len(n) > 1:
                                ins.append(n)
                        for insertions in ins:
                            nucs_here.discard(insertions)
                        final_code = get_non_iupac(frozenset(nucs_here))
                        if final_code is None:
                            print(frozenset(str(nucs_here)))
                            raise RuntimeError()
                        if ins:
                            ins.append(final_code)
                            yield yield_fstr.format('|'.join(ins))
                        else:
                            yield final_code
            else:
                yield self.seq[i].nuc

    @property
    def complexity(self):
        return sum(p.complexity for p in self.ambiguous.values())

    @property
    def ambiguous(self):
        return self._ambiguous

    @property
    def ref(self):
        return self._ref

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def seq(self):
        return self._seq

    @property
    def count(self):
        return self._count

    @count.setter
    def count(self, val):
        if not isinstance(val, (int, float)):
            raise RuntimeError("Count must be int or float! Got {} instead".format(val))

    @property
    def coverage(self):
        return int(float(self.count) / float(self.stop - self.start + 1))

    def __iter__(self):
        return iter(self._seq)

    def __len__(self):
        return len(self._seq)

    def __getitem__(self, index):
        return self._seq[index]


class ConsensusPosition(object):
    ambiguity_thresh = 1.0 / 3.0
    def __init__(self, pos, nuc, count, ref_call):
        self._pos = -1
        self._nuc = None
        self._count = -1
        self._ref_call = ''
        self._ambiguous = False
        self._analyze = True
        self._complexity = 0
        self.initialize(pos, nuc, count, ref_call)

    def initialize(self, pos, nuc, count, ref_call):
        self.pos = pos
        self.count = count
        self.ref_call = ref_call
        if isinstance(nuc, list):
            self.ambiguous = True
            self.nuc = list(map(str.upper, nuc))
        elif isinstance(nuc, str):
            self.nuc = nuc.upper()
        else:
            raise ValueError('Nucleotide position must be of type string')

    def flatten(self):
        if not self.ambiguous or not self.analyze:
            return self.ambiguous
        total = float(len(self.nuc))
        to_keep = set()
        counts = counter(self.nuc)
        for k, c in counts.items():
            if float(c) / total >= ConsensusPosition.ambiguity_thresh:
                to_keep.add(k)
        if len(to_keep) == 1:
            self.nuc = to_keep.pop()
            self.ambiguous = not self._nuc == self._ref_call
            self.count = counts[self.nuc]
            self._complexity = int(self.ambiguous)
        elif self.ref_call in to_keep:
            self.nuc = self.ref_call
            self.ambiguous = False
            self.count = counts[self.nuc]
        else:
            i = 0
            while i < len(self.nuc):
                if self.nuc[i] in to_keep:
                    i += 1
                else:
                    # remove the sequence
                    removed = self.nuc.pop(i)
                    # remove the count of nucleotides
                    self.count -= len(removed)
            bases = set(self.nuc)
            for n in bases:
                self._complexity += len(n)
        return self.ambiguous

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, val):
        if isinstance(val, int):
            self._pos = val
        else:
            raise ValueError('Nucleotide position should be of type int; got {} instead'.format(type(val)))

    @property
    def nuc(self):
        return self._nuc

    @nuc.setter
    def nuc(self, val):
        if isinstance(val, str):
            self._nuc = val
        elif isinstance(val, list):
            if not self.ambiguous:
                raise ValueError('You must specify that this position is an ambiguous position before setting ambiguous values')
            self._nuc = val
        else:
            raise ValueError('Nucleotide position must be string/list type; got {} instead'.format(type(val)))

    @property
    def count(self):
        return self._count

    @count.setter
    def count(self, val):
        if isinstance(val, int):
            self._count = val
        else:
            raise ValueError('Nucleotide count must be of type int; got {} instead'.format(type(val)))

    @property
    def ambiguous(self):
        return self._ambiguous

    @ambiguous.setter
    def ambiguous(self, val):
        if isinstance(val, bool):
            self._ambiguous = val
        else:
            raise ValueError('Ambiguity must be of type bool; got {} instead'.format(type(val)))

    @property
    def analyze(self):
        return self._analyze

    @analyze.setter
    def analyze(self, val):
        if isinstance(val, bool):
            self._analyze = val
        else:
            raise ValueError('Analyze attribute must be of type bool; got {} instead'.format(type(val)))

    @property
    def ref_call(self):
        return self._ref_call

    @ref_call.setter
    def ref_call(self, val):
        if isinstance(val, str):
            self._ref_call = val.upper()
        else:
            raise RuntimeError('Reference nucleotide must be string type got {} instead'.format(type(val)))

    @property
    def complexity(self):
        return self._complexity


def pileup_iterator(flname=None, flobj=None):
    if flobj:
        for line in _pileup_iterator(flobj):
            yield line
    elif flname:
        if not os.path.exists(flname):
            raise RuntimeError('Invalid pileup file')
        with open(flname, 'r') as f:
            for line in _pileup_iterator(f):
                yield line
    else:
        raise RuntimeError('No fileobj or filename provided for pileup scanning')

def _pileup_iterator(flobj):
    for line in flobj:
        split_line = line.strip().split('\t')
        while len(line) < 7:
            split_line.append('')
        yield split_line

def base_action(read_calls, index, final, match):
    return match.end()

def indel_action(read_calls, index, final, match):
    if read_calls[index] not in {"+", "-"}:
        return match.end()

    if read_calls[index] == '+':
        final[-1] += match.group(2).upper()
    # Add the number of deletions to the last sequence
    # We will need it to know how much to expand the
    # window by
    elif read_calls[index] == '-':
        final[-1] = (final[-1], int(match.group(1)))

    return match.start() + int(match.group(1)) + len("-" + match.group(1))

def _process_iter_re(read_calls, match_obj, index=0, final=None):
    if final is None:
        final = []
    upto = match_obj.start()
    while index < upto:
        final.append(read_calls[index].upper())
        index += 1
    return index, final

def process_iter_re(read_calls, matcher=None, match_objs=None, action=None):
    # the purpose of this function is to parse out the pile_up line
    # into discrete pieces:
    # E.g.
    #
    # ,  +1t  ,+1t  .+1T  .  .+1T  .+1T  .+1T  .  .+1T  .+1T
    # |    |             |  |
    # ^^^^^^             ^^^^
    #    |                  There are no insertions here
    # These two are associated
    #
    if action is None or not callable(action):
        action = base_action
    index = 0
    final = []
    if match_objs:
        for match in match_objs:
            index, final = _process_iter_re(read_calls, match, index=index, final=final)
            index = action(read_calls, index, final, match)
        while index < len(read_calls):
            final.append(read_calls[index].upper())
            index += 1
    elif matcher:
        for match in matcher.finditer(read_calls):
            index, final = _process_iter_re(read_calls, match, index=index, final=final)
            index = action(read_calls, index, final, match)
        while index < len(read_calls):
            final.append(read_calls[index].upper())
            index +=1
    return final

def sub_reference(string, reference_call):
    return string.replace('.', reference_call).replace(',', reference_call)

def sub_asterisk(string):
    return string.replace('*', '-')

def process_line(line):
    # All the information associated to a line
    ref = line[0]
    position = int(line[1])
    reference_call = line[2]
    read_count = int(line[3])
    read_calls = line[4]
    quality_scores = map(lambda x: ord(x)-33, line[5])
    mapping_quality = map(lambda x: ord(x)-33, line[6])
    ins = list(_insertions.finditer(read_calls))
    dels = list(_deletions.finditer(read_calls))
    if ins or dels:
        # Insertions are very straight forward,
        # Let's assume that the nucleotide position is a kmer of k > x
        # where x is the length of the insertion. The read count
        # is going to be the number of positions where there is an
        # insertion
        indels = list(ins)
        indels.extend(dels)
        indels.extend(m for m in _start_read.finditer(read_calls))
        indels.extend(m for m in _remove.finditer(read_calls))
        indels.sort(key=lambda x: x.start())
        extracted = process_iter_re(read_calls, match_objs=indels, action=indel_action)
        extracted_with_ref_calls = []
        if dels:
            for e in extracted:
                if isinstance(e, tuple):
                    extracted_with_ref_calls.append((e[0].replace('.', reference_call).replace(',', reference_call), e[1]))
                else:
                    extracted_with_ref_calls.append((e.replace(',', reference_call).replace('.', reference_call), 0))
        else:
            for e in extracted:
                new_e = sub_asterisk(sub_reference(e, reference_call))
                extracted_with_ref_calls.append(new_e)
        return bool(dels), ref, position, read_count, extracted_with_ref_calls, reference_call

    read_calls_list = []
    # Remove the dollar signs and other things
    # that don't tell us things we actually want to know
    read_calls = _remove.sub(r'', read_calls)
    read_calls_list = process_iter_re(read_calls, matcher=_start_read, action=base_action)
    # Replace all the dots and commas with the reference call
    for i in range(len(read_calls_list)):
        read_calls_list[i] = sub_reference(read_calls_list[i], reference_call)
        read_calls_list[i] = sub_asterisk(read_calls_list[i])
    if all(nuc == reference_call for nuc in read_calls_list):
        return False, ref, position, read_count, reference_call, reference_call
    else:
        return False, ref, position, read_count, read_calls_list, reference_call

def extend_detection_window(pile, processed_stack):
    for i, line in enumerate(processed_stack):
        if line[0]:
            max_jump = float('-inf')
            for deletion in line:
                if isinstance(deletion, tuple):
                    max_jump = max(max_jump, deletion[1])
            if len(processed_stack) - i - 1 < max_jump:
                for _ in range(len(processed_stack)-i-1):
                    processed_stack.append(process_line(next(pile)))

def detect_deletion_window(pile, ref, position, read_count, read_calls, current_consensus, ref_call):
    max_jump = max(read_calls, key=lambda x: x[1])[1]
    stack = []
    # In the weird off chance that the deletion is
    # at the end of the sequence but samtools  gives us a deletion
    # larger than the remainder of the data stream
    try:
        for _ in range(max_jump):
            stack.append(next(pile))
    except StopIteration:
        max_jump = len(stack)
    processed_stack = [process_line(line) for line in stack]
    # Check to see if we have any extra deletions in the the window we just pulled
    extend_detection_window(pile, processed_stack)
    all_lines = [read_calls]
    all_lines.extend(line[4] for line in processed_stack)
    all_lines_fixed = []
    for i, l in enumerate(all_lines):
        if isinstance(l, list):
            all_lines_fixed.append(list())
            for n in l:
                if isinstance(n, tuple):
                    all_lines_fixed[-1].append(n[0])
                else:
                    all_lines_fixed[-1].append(n)
        else:
            if i == 0:
                count = read_count
            else:
                count = processed_stack[i-1][3]
            all_lines_fixed.append(repeat(l, count))
    #
    # all_lines represents a window: [
    #   POS 108 -> [nuc, nuc, ... , nuc, nuc, nuc]
    #   POS 109 -> [nuc, nuc, ... , nuc, nuc]
    #   POS 110 -> [nuc, nuc, ... , nuc]
    #   POS 111 -> [nuc, nuc, ... , nuc, nuc]
    # ]
    #
    #
    # The first list in all_lines is the read_calls this function was called with
    # so the first position prior to the deletion window
    # and the last list being the last position in the deletion window
    #
    sub_seq = []
    for window in zip_longest(*all_lines_fixed, fillvalue='-'):
        # If bowtie says that this is a deletion window
        # Then it better be consistent
        # Namely, this happens:
        #   POS 108 -> .-3TTA
        #   POS 109 -> A
        #   POS 110 -> *
        #   POS 111 -> T
        #
        # Clearly this is stupid, mpileup says there's a deletion here
        # but that is not reflected in subsequent reads
        #
        # The opposite will also happen where:
        #
        #   POS 108 ->  .
        #   POS 109 ->  *
        #   POS 110 ->  *
        #   POS 111 ->  *
        #
        # Clearly a deletion event of 3 was not flagged
        # and yet all reads afterwards are deleted
        if all(x=='-' for x in window[1:]):
            sub_seq.append(tuple(window))
        elif all(x!='-' for x in window[1:]):
            sub_seq.append(tuple(window))
    # We want to get the counts of the deletion window
    # meaning how many reads support a certain substring
    counts = counter(sub_seq)
    to_remove = []
    for seq, count in counts.items():
        if float(count) / len(sub_seq) < ConsensusPosition.ambiguity_thresh:
            to_remove.append(seq)
    for rm in to_remove:
        del counts[rm]
    # If we reduce the complexity down
    # to just one substring, hooray!
    if len(counts) == 0:
        first_cp = ConsensusPosition(
            position,
            '-',
            0,
            ref_call
        )
        first_cp.ambiguous = True
        first_cp.analyze = False
        current_consensus.add_nuc(
            ref,
            first_cp.pos,
            first_cp,
            first_cp.count,
            ref_call
        )
        for i in range(len(processed_stack)):
            cp = ConsensusPosition(
                processed_stack[i][2],
                '-',
                0,
                processed_stack[i][5]
            )
            cp.ambiguous = True
            cp.analyze = False
            current_consensus.add_nuc(
                processed_stack[i][1],
                cp.pos,
                cp,
                cp.count,
                processed_stack[i][5]
            )
        return 0
    elif len(counts) == 1:
        final_seqs = list(counts.keys())[0]
        count = list(counts.values())[0]
        current_consensus.add_nuc(ref, position, final_seqs[0], count, ref_call)
        for i, next_nuc in enumerate(final_seqs[1:]):
            current_consensus.add_nuc(
                ref,
                processed_stack[i][2],
                next_nuc,
                0 if next_nuc=='-' else count,
                processed_stack[i][5]
        )
        return 0
    else:
        i = len(sub_seq)-1
        while i >= 0:
            if sub_seq[i] not in counts:
                sub_seq.pop(i)
            i-=1
        final_seqs = list(map(list, zip(*sub_seq)))
        first_cp = ConsensusPosition(
            position,
            final_seqs[0],
            sum(x!='-' for x in final_seqs[0]),
            ref_call
        )
        current_consensus.add_nuc(
            ref,
            first_cp.pos,
            first_cp,
            first_cp.count,
            ref_call
        )
        for i, next_nuc in enumerate(final_seqs[1:]):
            cp = ConsensusPosition(
                processed_stack[i][2],
                next_nuc,
                sum(x!='-' for x in next_nuc),
                processed_stack[i][5]
            )
            current_consensus.add_nuc(
                processed_stack[i][1],
                cp.pos,
                cp,
                cp.count,
                processed_stack[i][5]
            )
        return 0

def build_consensus(pileup_file, ambiguity_threshold):
    reference = None
    current_consensus = ConsensusSequence()
    position_offset = 0
    pile = pileup_iterator(flname=pileup_file)
    ConsensusPosition.ambiguity_thresh = ambiguity_threshold
    while pile:
        try:
            next_line = next(pile)
        except StopIteration:
            return
        dels, ref, position, read_count, read_calls, reference_call = process_line(next_line)
        if position < current_consensus.stop:
            # We've started on the next reference in the pileup file
            yield reference, current_consensus
            current_consensus = ConsensusSequence()
            position_offset = 0
        reference = ref
        if not read_count:
            current_consensus.add_nuc(reference, position + position_offset, '-', read_count, reference_call)
            continue
        if dels:
            position_offset += detect_deletion_window(
                pile,
                ref, position,
                read_count,
                read_calls,
                current_consensus,
                reference_call
            )
            continue
        current_consensus.add_nuc(reference, position + position_offset, read_calls, read_count, reference_call)

def build_sequences(pileup_file, ambiguity_threshold, min_coverage, max_complexity):
    final_sequences = {}
    final = []
    for reference, c_sequence in build_consensus(pileup_file, ambiguity_threshold):
        if c_sequence.coverage < min_coverage:
            continue
        if 'stx2a' in reference:
            pause = 1
        c_sequence.flatten()
        new_ref = reference.split('|')[0].split('_')[0]
        if new_ref in final_sequences:
            final_sequences[new_ref].append(c_sequence)
        else:
            final_sequences[new_ref] = [c_sequence]
    for ref, seqs in final_sequences.items():
        if not seqs:
            continue
        seqs.sort(key=lambda x: x.complexity)
        if seqs[0].complexity <= max_complexity:
            final.append(seqs[0])
    return final


class SamtoolsFilter(object):
    """
    Use this class to create filter objects for sam files to pass to
    sam view.
    These are the default alignment flags:
    0x1     PAIRED  paired-end (or multiple-segment) sequencing technology
    0x2     PROPER_PAIR each segment properly aligned according to the aligner
    0x4     UNMAP   segment unmapped
    0x8     MUNMAP  next segment in the template unmapped
    0x10    REVERSE SEQ is reverse complemented
    0x20    MREVERSE    SEQ of the next segment in the template is reverse complemented
    0x40    READ1   the first segment in the template
    0x80    READ2   the last segment in the template
    0x100   SECONDARY   secondary alignment
    0x200   QCFAIL  not passing quality controls
    0x400   DUP PCR or optical duplicate
    0x800   SUPPLEMENTARY   supplementary alignment
    """
    flags = {
        'paired': 0x1,
        'proper_pair': 0x2,
        'unmapped': 0x4,
        'mate_unmapped': 0x8,
        'reversed': 0x10,
        'mate_reversed': 0x20,
        'read_1': 0x40,
        'read_2': 0x80,
        'secondary': 0x100,
        'qc_fail': 0x200,
        'duplicate': 0x400,
        'supplementary': 0x800
    }
    def __init__(self, filter_flags=None, require_flags=None):
        self._require_flag = 0x0
        self._filter_flag = 0x0
        if filter_flags is not None:
            if isinstance(filter_flags, (set, frozenset)):
                filter_flags = list(filter_flags)
            elif not isinstance(filter_flags, list):
                raise TypeError('Filter flags must be list type!')
            self.filter(*filter_flags)

        if require_flags is not None:
            if isinstance(require_flags, (set, frozenset)):
                require_flags = list(require_flags)
            elif not isinstance(require_flags, list):
                raise TypeError('Require flags must be list type!')
            self.require(*require_flags)

    def require(self, *args):
        for flag in set(args):
            if flag not in SamtoolsFilter.flags:
                raise RuntimeError('Inappropriate flag provided: {}'.format(flag))
            self._require_flag |= SamtoolsFilter.flags[flag]

    def filter(self, *args):
        for flag in set(args):
            if flag not in SamtoolsFilter.flags:
                raise RuntimeError('Inappropriate flag provided: {}'.format(flag))
            self._filter_flag |= SamtoolsFilter.flags[flag]

    def reset(self, require=False, flter=False):
        if require:
            self._require_flag = 0x0
        if flter:
            self._filter_flag = 0x0

    @property
    def require_flag(self):
        return str(self._require_flag)

    @property
    def filter_flag(self):
        return str(self._filter_flag)


def valid_dir(path):
    # Check to make sure the directory exists
    if os.path.exists(path):
        return

    if not os.path.isdir(path):
        os.makedirs(path)

def check_dir(path):
    return os.path.isdir(path)

def full_path(path):
    return os.path.realpath(
                os.path.expanduser(
                    os.path.expandvars(path)
                ))

class EcStxTyper(Algorithm[EcStxTyperInputs, EcStxTyperOutputs]):

    outputs_class = EcStxTyperOutputs

    def execute_stub(self):
        working_dir = Path("tmp")
        working_dir.mkdir(exist_ok=True)
        cmd = [f"{CONDA_BIN_DIR}/bowtie2"]
        process_stdout, _ = run_external(cmd, self.logger, True, working_dir)
        print(process_stdout)

    def execute_implementation(self):
        ref_file = self.fmt_reference_file()
        ref_index = self.run_bowtie_index(ref_file)
        sam_file = self.run_bowtie(ref_index)
        samview_file = self.run_samview(sam_file)
        sorted_bam_file = self.run_bamsort(samview_file)
        pileup_file = self.run_pileup(ref_file, sorted_bam_file)
        self.ecstxtyper_result_out = self.generate_results(pileup_file)
        self.outputs.ecstxtyper_result_out.to_file()
        self.set_results()

    def fmt_reference_file(self):
        STXTarget = namedtuple("STXTarget", ["locus", "allele", "accession", "sequence"])
        def sequence_parser(header, sequence, sep="_"):
            info = header.split(sep)
            while len(info) < 4:
                info.append("")
            return STXTarget(locus=info[0] + info[1],
                             allele=info[2],
                             accession=info[3],
                             sequence=sequence
                             )

        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        db_path = self.inputs.holotoxins
        seq_db = DbInfo(db_path, seq_parser=sequence_parser)
        outfile = working_dir / "stxrefs.fasta"
        seq_db.export_sequences(outfile)
        return working_dir / "stxrefs.fasta"

    def run_bowtie_index(self, ref_file: Path, name=""):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        if name:
            ref_name = name
        else:
            ref_name = os.path.basename(ref_file).split('.')[0]

        idx_dir = full_path(os.path.join('bowtie', ref_name, 'index'))
        valid_dir(idx_dir)
        cmd = ["bowtie2-build"]
        cmd += [ref_file]
        cmd += [idx_dir]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return idx_dir

    def run_bowtie(self, index_dir: Path):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        cmd = ["bowtie2"]
        cmd += ["-x", index_dir]
        cmd += ["-p", str("1")]
        cmd += ["--reorder"]
        cmd += ["--local"]
        cmd += ["--sensitive-local"]
        cmd += ["--no-unal"]
        cmd += ["--all"]
        cmd += ["-1", str(self.inputs.read1)]
        cmd += ["-2", str(self.inputs.read2)]
        cmd += ["-S", "output.sam"]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return working_dir / "output.sam"

    def run_samview(self, sam_file: Path):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        cmd = ["samtools"]
        cmd += ["view"]
        cmd += ["-b"]
        cmd += ["-S"]
        cmd += ["-o", "SVoutput.bam"]
        cmd += [str(sam_file)]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return working_dir / "SVoutput.bam"

    def run_bamsort(self, samview_file: Path):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        bam_dir, bam_name = os.path.split(samview_file)
        sorted_name = bam_name.split('.')[0] + '.sorted'
        bam_sorted_path = full_path(os.path.join(bam_dir, sorted_name))
        intermediate_name = bam_name.split('.')[0] + '_intermediate'
        cmd = ["samtools"]
        cmd += ["sort"]
        cmd += ["-T", intermediate_name]
        cmd += ["-o", bam_sorted_path]
        cmd += [str(samview_file)]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return bam_sorted_path

    def run_pileup(self, ref_file: Path, sorted_bam: Path):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        pileup_name = "stx.pup"
        default_filter = SamtoolsFilter()
        default_filter.filter(
            "unmapped",
            "qc_fail",
            "duplicate"
        )
        cmd = ["samtools", "mpileup"]
        cmd += ["-aa"]
        cmd += ["-f", ref_file]
        cmd += ["--ff", default_filter.filter_flag]
        cmd += ["-s", sorted_bam]
        cmd += ["-o", str(pileup_name)]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return working_dir / str(pileup_name)

    def generate_results(self, pileup_file: Path):

        def results_parser(dbinfo, results):
            results_out = {
                'results': {},
                'extra': []
            }
            for info in dbinfo.sequences.values():
                gene_name = info.locus
                results_out['results'][gene_name] = False

            for seq in results:
                reference_info = dbinfo.sequences.get(seq.ref.split('|')[0], None)
                results_out['results'][reference_info.locus] = True

                results_out['extra'].append({
                    'locus': reference_info.locus,
                    'allele': reference_info.allele,
                    'accession': reference_info.accession,
                    'sequence': ''.join(seq.get_fragment()),
                    'coverage': str(seq.coverage),
                    'complexity': str(seq.complexity)
                })
            return results_out

        STXTarget = namedtuple("STXTarget", ["locus", "allele", "accession", "sequence"])
        def sequence_parser(header, sequence, sep="_"):
            info = header.split(sep)
            while len(info) < 4:
                info.append("")
            return STXTarget(locus=info[0] + info[1],
                             allele=info[2],
                             accession=info[3],
                             sequence=sequence
                             )

        db_path = self.inputs.holotoxins
        seq_db = DbInfo(db_path, seq_parser=sequence_parser)

        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)

        found_sequences = build_sequences(pileup_file, float(0.333333), float(10.0), float(3.0))
        results_out = seq_db.results_parser(found_sequences, f=results_parser)

        self.outputs.ecstxtyper_result_out.content = results_out

    def set_results(self):
        self.result = {
            "ecstxtyper_result_out": self.inputs.publish_dir + str(self.outputs.ecstxtyper_result_out.path)
        }