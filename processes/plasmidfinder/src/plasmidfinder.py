from collections import defaultdict, Counter, namedtuple
from csv import DictReader
from os import getenv, makedirs, listdir
from os.path import (
    exists,
    join,
    isdir,
    splitext,
    dirname
)
from pathlib import Path
from shutil import move, rmtree
from itertools import combinations

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.base.inputs import OrganismInput
from ngs_pipeline_lib.tools.tools import gunzip_file

from src.inputs import PlasmidFinderInputs
from src.outputs import PlasmidFinderOutputs
from src.config import PLASMIDFINDER_SETTINGS


class GenotypeHit(object):

    def __init__(self):
        self.reference_id = ''
        self.reference_len = ''
        self.query_id = ''
        self._query_start = 0
        self._query_stop = 0
        self._reference_start = 0
        self._reference_stop = 0
        self.forward = True
        self.identity = 0.
        self.absolute_len = 0
        self.relative_len = 0
        self.query_len = 0
        self.num_mismatches = 0
        self.num_gap_opens = 0
        self.query_seq = None
        self.reference_seq = None
        self.bitscore = 0
        self.evalue = 0
        self.full_match = False

    # ************* IMPORTANT *************
    # 
    # The following methods are used to cope
    # with using 1 indexed BLAST indices.
    # If you are going to use the setters
    # make sure you set it with the 0 indexed
    # value, which is reasonable since you are getting
    # the 0 indexed values from the below
    # property methods to begin with.
    # We did this so I can have good feels about this
    # datastructure (thanks).
    # 

    @property
    def query_start(self):
        return self._query_start - 1

    @query_start.setter
    def query_start(self, val):

        self._query_start = val + 1

    @property
    def query_stop(self):
        return self._query_stop - 1

    @query_stop.setter
    def query_stop(self, val):
        self._query_stop = val + 1

    @property
    def reference_start(self):
        return self._reference_start - 1

    @reference_start.setter
    def reference_start(self, val):
        self._reference_start = val + 1

    @property
    def reference_stop(self):
        return self._reference_stop - 1

    @reference_stop.setter
    def reference_stop(self, val):
        self._reference_stop = val + 1

    @staticmethod
    def from_blast(line):
        # query_id  ref_id|len  iden    alignmentlen    mismatches
        #   gapopens    qstart  qstop   refstart    refstop evalue
        #   bitscore    qlen    queryseq    refseq 

        # Create the object
        hit = GenotypeHit()

        # Split the line
        parts = line.split('\t')

        # Parse the results
        hit.query_id = parts[0]
        hit.reference_id = parts[1]
        hit.reference_len = int((parts[1]).split("|")[1])
        hit.identity = float(parts[2])
        hit.absolute_len = int(parts[3])
        hit.num_mismatches = int(parts[4])
        hit.num_gap_opens = int(parts[5])

        # NOTE:
        # The alignments come in 1-indexed
        hit._query_start = int(parts[6])
        hit._query_stop = int(parts[7])
        hit._reference_start = int(parts[8])
        hit._reference_stop = int(parts[9])

        hit.evalue = float(parts[10])
        hit.bitscore = float(parts[11])

        hit.query_len = int(parts[12])
        hit.coverage = float((hit._query_stop-hit._query_start+1)/(hit.reference_len))
        # Not all of the genotyping needs to 
        # have the sequences returned
        if len(parts) == 15:
            hit.query_seq = parts[13]
            hit.reference_seq = parts[14]

        # Check to see if the sequence has been reversed
        # If for whatever reason the hit is one nucleotide
        # long, this will be false
        hit.forward = hit._reference_start < hit._reference_stop

        if not hit.forward:
            hit._reference_start, hit._reference_stop = \
                hit._reference_stop, hit._reference_start

        # Change the percentages into decimals,
        # because why not
        hit.identity = hit.identity / 100.

        # Make sure to avoid a divide by zero error
        if hit.query_len:
            hit.relative_len = float(hit.absolute_len) / float(hit.query_len)

        if not hit.num_gap_opens and hit.identity == 1.0 and \
            hit.relative_len == 1.0:
            hit.full_match = True

        return hit
    
class GenotypeResults(object):

    def __init__(self, filename):
        self._hits = []
        self.load_hits(filename)

    def load_hits(self, filename):
        if isinstance(filename, str):
            if exists(filename):
                with open(filename, 'r') as f:
                    for line in f.readlines():
                        self._hits.append(GenotypeHit.from_blast(line))

            else:

                raise RuntimeError('Provide file path is not a real'
                    ' path for alignment results parsing')

        else:
            # We were given a file handle
            for line in filename.readlines():

                self._hits.append(GenotypeHit.from_blast(line))

    @property
    def hits(self):
        return self._hits


SequenceInfo = namedtuple('SequenceInfo', ['locus', 'allele', 'accession', 'sequence', 'other'])
LocusInfo = namedtuple('LocusInfo', ['locus', 'note', 'antibiotic', 'other'])

class PlasmidFinder(Algorithm[PlasmidFinderInputs, PlasmidFinderOutputs]):

    outputs_class = PlasmidFinderOutputs

    def execute_stub(self):
        pass

    def execute_implementation(self):
        self.working_dir = Path("./")
        blast_out = Path(f"{self.working_dir}_{self.inputs.sample_id}_plasmidfinder_blast.txt")
        genus = self.inputs.organism.genus
        settings = PLASMIDFINDER_SETTINGS[genus]
        
        # make blastdb each time this job is called
        sequences, concat_file_path = self.concatenate_fasta(dirpath=self.inputs.reference)
        blast_db_path = "db_concat"
        self.make_blast_db(
            concat_file_path,
            blast_db_path
        )

        # Read assembly file
        assembly = gunzip_file(self.inputs.assembly)
        assembly_dict: dict = {}
        with open(assembly) as reader:
            for header, seq in SimpleFastaParser(reader):
                assembly_dict[header] = seq

        # run blastn
        results = self.run_blast(
            query = assembly,
            db_path = blast_db_path,
            percent_identity = settings["PF_PERCENT_IDENTITY"],
            outfile = blast_out
        )

        #
        results = self.results_parser(
            results,
            sequences,
            blast_coverage=settings["PF_MIN_RELATIVE_COVERAGE"],
            min_overlap = settings["PF_MIN_MERGE_OVERLAP"]
        )

        bheader = "Query\tSubject\tPercent Identity\tAlignment Length\
            \tMismatches\tGap Opens\tQuery Start\tQuery End\tSubject Start\
            \tSubject End\tEvalue\tBit Score\tQuery Length\n"
        with open(blast_out, 'r') as bfile: 
            bcontent = bheader + bfile.read()

        self.outputs.plasmidfinder.content = results
        self.outputs.plasmidfinder.to_file()
        self.outputs.blastout.content = bcontent
        

    def fasta_iterator(self, file_name):
        with open(file_name, 'r') as f:
            sequence_parts = []
            key = ''
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line[0] == '>':
                    if key:
                        full_sequence = ''.join(sequence_parts).upper()
                        yield (key, full_sequence)
                    key = line[1:].split()[0]
                    sequence_parts = []
                else:
                    sequence_parts.append(line)
            if key:
                full_sequence = ''.join(sequence_parts).upper()
                yield (key, full_sequence)

    def is_fasta(self, path: str) -> bool:
        _FASTAEXTS = ['.fna', '.fasta', '.fsa']
        if not exists(path):
            return False
        file_extension = splitext(path)[1]
        return file_extension in _FASTAEXTS


    def parse_fasta(self, file_name, rename=False):
        if self.is_fasta(file_name):
            fasta_sequences = {}
            if rename:
                for i, (name, sequence) in enumerate(self.fasta_iterator(file_name), 1):
                    new_name = 'contig_' + str(i)
                    fasta_sequences[new_name] = sequence
                return fasta_sequences
            else:
                fasta_sequences.update(self.fasta_iterator(file_name))
                return fasta_sequences
    
    def sequence_parser(self, header, sequence, sep=':'):
        parts = header.split(sep)
        parts = list(map(str.strip, parts))
        while len(parts) < 4:
            parts.append('')

        return SequenceInfo(locus = parts[0],allele = parts[1],accession = parts[2],sequence = sequence,other = parts[3])

    def valid_dir(self, path):
        # Check to make sure the directory exists
        if exists(path):
            return
        if not isdir(path):
            makedirs(path)

    def export_sequences(self, file_path: Path, seq_dict: dict) -> Path:
        self.valid_dir(dirname(file_path))
        with open(file_path, 'w') as f:
            for seq_id, seq_info in seq_dict.items():
                # f.write(f">{seq_id}\n{seq_info}\n")
                ostr = '>{}|{}\n{}\n'
                f.write(ostr.format(seq_id,len(seq_info),seq_info))

    def concatenate_fasta(self, dirpath: str|Path) -> tuple[dict[str,str],Path]:
        """
        :param dirpath: finderdb directory, may contain multiple fasta files
        :return: sequences dictionary
        """
        dirpath=Path(dirpath)
        sequences = {}
        sequence_counts = defaultdict(dict)
        allele_id_template = '{}_{}'
        allele_id_template_i = '{}_{}_{}'
        separator = '-'

        for seq_file in listdir(dirpath):
            # skip previous run's concatenated sequences
            if seq_file == "concat_sequences.fasta":
                continue
            file_path = join(dirpath, seq_file)
            if not self.is_fasta(file_path):
                continue
            current_sequences = self.parse_fasta(file_path)
            for seq_id, sequence in list(current_sequences.items()):
                seq_info = self.sequence_parser(seq_id, sequence)
                if seq_info is None:
                    continue
                allele_id = allele_id_template.format(
                    seq_info.locus, seq_info.allele
                )
                if allele_id in sequence_counts:
                    new_id = allele_id_template_i.format(
                        seq_info.locus,
                        seq_info.allele,
                        len(sequence_counts[allele_id])
                    )
                    sequence_counts[allele_id][new_id] = True
                    sequences[new_id] = seq_info.sequence
                else:
                    sequence_counts[allele_id][allele_id] = True
                    sequences[allele_id] = seq_info.sequence
        # Write the sequence dictionary to file
        concat_file = join(dirpath, 'concat_sequences.fasta')

        if not exists(concat_file):
            self.export_sequences(file_path=concat_file, seq_dict=sequences)
        return sequences, concat_file

    def make_blast_db(
        self,
        concat_file: Path,
        db_out_path: Path
    ) -> Path:
        cmd = ["mamba run"]
        cmd = ["makeblastdb"]
        cmd += ["-in", concat_file]
        cmd += ["-dbtype", "nucl"]
        cmd += ["-out", db_out_path]
        
        blast_db_process_stdout, _ = run_external(
            cmd,
            self.logger,
            use_mamba_env=False,
            text=True,
            working_dir=self.working_dir
            )
        self.logger.info(blast_db_process_stdout)

    def run_blast(self, query: Path, db_path: Path, percent_identity: float, outfile: Path) -> GenotypeResults:
        cmd = ["mamba run"]
        cmd = ["blastn"]
        cmd += ["-query", str(query)]
        cmd += ["-db", db_path]
        cmd += ["-task", "dc-megablast"]
        cmd += ["-out", str(outfile)]
        cmd += ["-perc_identity", str(100.0*percent_identity)]
        cmd += ["-outfmt", "6 std qlen"]
        cmd += ["-max_target_seqs",str(1000000)]
        cmd += ["-dust", "no"]
        blast_process_stdout, _ = run_external(
            cmd,
            self.logger,
            use_mamba_env=False,
            text=True,
            working_dir=self.working_dir
        )
        self.logger.info(blast_process_stdout)
        return GenotypeResults(str(outfile))
    


    
    def results_parser(self, results:GenotypeResults, sequences:dict, blast_coverage: float, min_overlap: float ) -> dict:
        
        results_out = {
            'results': {},
            'extra': []
        }
        ref_id_dict = {}

        #get best hit with same ref id based on coverage*identity
        for hit in results.hits: 
            if hit.reference_id in ref_id_dict.keys():
                if (hit.coverage*hit.identity) > ((ref_id_dict[hit.reference_id]).coverage)*((ref_id_dict[hit.reference_id]).identity):
                    del ref_id_dict[hit.reference_id]
                    ref_id_dict[hit.reference_id] = hit
            else:
                ref_id_dict[hit.reference_id] = hit

        #get best hit with different ref ids based on identity
        to_remove = set()
        for hit1 in ref_id_dict.values():
            for hit2 in ref_id_dict.values():
                if hit1.coverage >= blast_coverage and hit2.coverage >= blast_coverage: 
                    if hit1.query_id == hit2.query_id and hit1 != hit2:
                        overlap = min(hit1.query_stop, hit2.query_stop) - \
                        max(hit1.query_start, hit2.query_start) + 1
                        
                        overlap = max(0, overlap)

                        hit1_length = hit1.query_stop - hit1.query_start + 1
                        hit2_length = hit2.query_stop - hit2.query_start + 1
                        
                        if overlap >= min_overlap * min(hit1_length, hit2_length):
                            if hit1.identity > hit2.identity:
                                to_remove.add(hit2.reference_id)
                            if hit1.identity < hit2.identity:
                                to_remove.add(hit1.reference_id)
        
        for rm in to_remove:
            del ref_id_dict[rm]

        for hit in ref_id_dict.values():
            locus = hit.reference_id.split("_")[0]
            allele = hit.reference_id.split("_")[1]
            gene_name = locus
            if hit.coverage >= blast_coverage:   
                results_out['results'][gene_name] = 1 ## 
                hit_information =[
                    {
                        'locus': locus,
                        'identity': hit.identity,
                        'allele': allele,
                        'hits': [
                            {
                                'contig_id': hit.query_id,
                                'reference_id':hit.reference_id,
                                'query_start': hit.query_start,
                                'query_stop': hit.query_stop,
                                'reference_start': hit.reference_start,
                                'reference_stop': hit.reference_stop,
                                'full_match': hit.full_match
                            }
                        ]   
                    }
                ]
                results_out['extra'].extend(hit_information)
            else:
                if gene_name in results_out['results']:
                    continue
                results_out['results'][gene_name] = False
        for sequence_info in sequences.items():
            locus = sequence_info[0].split("_")[0]
            allele = sequence_info[0].split("_")[1]
            gene_name = locus
            if gene_name not in results_out['results']:
                results_out['results'][gene_name] = 2 ##

        
        return results_out
