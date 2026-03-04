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
import copy

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.base.inputs import OrganismInput
from ngs_pipeline_lib.tools.tools import gunzip_file

from src.inputs import PathotypeFinderInputs
from src.outputs import PathotypeFinderOutputs
from src.config import SETTINGS

_PATHOTYPES = {
    'STEC': "Shiga toxin-producing Escherichia coli",
    'ETEC': "Enterotoxigenic Escherichia coli",
    'EPEC': "Enteropathogenic Escherichia coli",
    'EPEC (typical)': "Enteropathogenic Escherichia coli (typical)",
    'EAEC': "Enteroaggregative Escherichia coli",
    'EIEC/Shigella': "Enteroinvasive Escherichia coli",
    'DAEC': "Diffusely adherent Escherichia coli",
    'STEC/EAEC': "Hybrid shiga toxin-producing / enteroaggregative Escherichia coli",
}

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
GenotypeRegion = namedtuple('GenotypeRegion', ['coverage', 'identity', 'locations'])

class PathotypeFinder(Algorithm[PathotypeFinderInputs, PathotypeFinderOutputs]):

    outputs_class = PathotypeFinderOutputs

    def execute_stub(self):
        pass

    def execute_implementation(self):
        self.working_dir = Path("./")
        blast_out = Path(f"{self.working_dir}_{self.inputs.sample_id}_pathotypefinder_blast.txt")
        genus = self.inputs.organism.genus
        settings = SETTINGS[genus]

        # make blastdb each time this job is called
        sequences, concat_file_path = self.concatenate_fasta(dirpath=self.inputs.reference)
        blast_db_path = "db_concat"
        self.make_blast_db(
            concat_file_path,
            blast_db_path
        )

        # Read assembly file
        assembly = gunzip_file(self.inputs.assembly)

        # run blastn
        results = self.run_blast(
            query = assembly,
            db_path = blast_db_path,
            percent_identity = settings["PERCENT_IDENTITY"],
            outfile = blast_out
        )

        #parse results
        results_out = self.results_parser(
            results,
            sequences,
            blast_coverage = settings["MIN_RELATIVE_COVERAGE"],
            min_overlap = settings["MIN_MERGE_OVERLAP"],
            percent_identity = settings["PERCENT_IDENTITY"],
            search_fragments = settings["SEARCH_FRAGMENTS"]
        )

        final_pathotypes = {
            'results' : set(),
            'extras' : []
        }

        present_genes = [gene for gene, boolean in results_out['results'].items() if boolean]

        final_pathotypes['results'].update(self.get_pathotypes(present_genes))

        for typ in final_pathotypes['results']:
            final_pathotypes['extras'].extend(_PATHOTYPES[typ] for typ in final_pathotypes['results'])

        final_pathotypes['results'] = list(final_pathotypes['results'])


        bheader = "Query\tSubject\tPercent Identity\tAlignment Length\
            \tMismatches\tGap Opens\tQuery Start\tQuery End\tSubject Start\
            \tSubject End\tEvalue\tBit Score\tQuery Length\n"
        with open(blast_out, 'r') as bfile: 
            bcontent = bheader + bfile.read()

        self.outputs.pathotypefinder_geno.content = results_out
        self.outputs.pathotypefinder_geno.to_file()
        self.outputs.blastout.content = bcontent

        self.result = final_pathotypes
    
    def concatenate_fasta(self, dirpath: str|Path) -> tuple[dict[str,str],Path]:
        """
        :param dirpath: may contain multiple fasta files
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

    #Needed for concatenate_fasta(), parse_fasta()
    def is_fasta(self, path: str) -> bool:
        _FASTAEXTS = ['.fna', '.fasta', '.fsa']
        if not exists(path):
            return False
        file_extension = splitext(path)[1]
        return file_extension in _FASTAEXTS

    #Needed for concatenate_fasta()
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

    #Needed for concatenate_fasta(), parse_fasta()
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

    #Needed for concatenate_fasta()
    # **pathotype db specific sequence parser
    def sequence_parser(self, header, sequence, sep = ':'):
        # etec.fasta has has a pipe in the references file
        # that's why we have to do this
        parts = header.split('|')[0].split(sep)
        while len(parts) < 4:
            parts.append('')

        return SequenceInfo(
            locus = parts[0],
            allele = parts[1],
            accession = parts[2],
            sequence = sequence,
            other = ''
        )

    #Needed for concatenate_fasta()
    def export_sequences(self, file_path: Path, seq_dict: dict) -> Path:
        self.valid_dir(dirname(file_path))
        with open(file_path, 'w') as f:
            for seq_id, seq_info in seq_dict.items():
                ostr = '>{}|{}\n{}\n'
                f.write(ostr.format(seq_id,len(seq_info),seq_info))

    #Needed for concatenate_fasta(), export_sequences()
    def valid_dir(self, path):
        # Check to make sure the directory exists
        if exists(path):
            return
        if not isdir(path):
            makedirs(path)

    def make_blast_db(
        self,
        concat_file: Path,
        db_out_path: Path
    ) -> Path:
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

    # Parse results. Search fragments implemented
    def results_parser(self, results:GenotypeResults, sequences:dict, blast_coverage: float, min_overlap: float, percent_identity: float, search_fragments: bool ) -> dict:
        
        #Check if hit is at edge when search fragments is true
        def at_edge(hit):

            return hit.query_start <= 25 or hit.query_len-hit.query_stop <= 25
             
        results_out = {
            'results': {},
            'extra': []
        }

        # Organize by reference
        hits_by_ref = {}
        for hit in results.hits:
            reference = hit.reference_id
            if reference not in hits_by_ref:
                hits_by_ref[reference] = []
            hits_by_ref[reference].append(hit)

        # Iterate through each ref's hits to validate and check for fragments, if applicable
        best_hits = {}
        for reference, hits, in hits_by_ref.items():
            to_check = []
            valid_hits = []
            for hit in hits:
                if hit.coverage >= blast_coverage and hit.identity >= percent_identity:
                    geno_region = GenotypeRegion(
                        coverage = hit.coverage,
                        identity = hit.identity,
                        locations = [hit]
                    )
                    valid_hits.append(geno_region)

                else:
                    to_check.append(hit)

            #Search for fragments
            if search_fragments:
                
                # Get a list of all of the edge hits
                to_check2 = list(filter(at_edge, to_check))

                for hits in combinations(to_check2, 2):
                    # Create a mask of the reference to see what
                    # amount of the reference is covered

                    hit1, hit2 = hits
                    mask = [0.] * hit1.reference_len
                    
                    # For each of the hits set the values
                    for hit in hits:
                        # Create the slice
                        # Plus 1 since list slicing is exclusive e.g. [1, 10) 
                        mask_slice = slice(
                            hit.reference_start,
                            hit.reference_stop+1,
                            None
                        )

                        # Size of the slice
                        size = hit.reference_stop - hit.reference_start + 1

                        # Set the mask slice
                        mask[mask_slice] = [hit.identity] * size

                    # Get the coverage of the two fragments
                    filtered_mask = list(filter(None, mask))
                    coverage = float(len(filtered_mask)) / float(hit1.reference_len)

                    # Get the real coverage, incase calculated
                    # coverage is greater than 1
                    coverage = min(1.0, coverage)

                    # If the coverage of the paired hits is good
                    # then check the identity
                    if coverage >= blast_coverage:

                        # The identity of the fragments is the best
                        # identities of both of the fragments
                        # divided by its relative coverage
                        # which comes out to be just the number of
                        # items in the mask greater than 0
                        identity = sum(mask) / (coverage * float(hit1.reference_len))

                        if identity >= percent_identity:
                            geno_region = GenotypeRegion(
                                coverage = coverage,
                                identity = identity,
                                locations=[hit1,hit2]
                            )
                            valid_hits.append(geno_region)

            # If there are no valid hits, go to next reference
            # Else, keep best hit for this reference
            if not len(valid_hits):
                continue
            else:
                valid_hits.sort(key=lambda x: float(x.coverage)*float(x.identity))
                best_hits[reference] = valid_hits.pop()

        # Find overlapping hits and fragments
        to_remove = set()
        checked_regions = [] #add to list to not compare twice
        for region1 in best_hits.values():
            checked_regions.append(region1)
            for region2 in best_hits.values():

                #skip if region has already been compared in the first loop
                if region2 in checked_regions:
                    continue

                if region1.coverage >= blast_coverage and region2.coverage >= blast_coverage: 
                    
                    # Find the regions that should be checked for overlap and put in pairs
                    # With fragments, there can be mulitple hits per reference so need to pair up
                    # the different fragments based on query id

                    compare = [] #will hold pairs to compare and check for overlap
                    for hit1 in region1.locations:
                        
                        hit1.identity = region1.identity
                        combine = []
                        combine.append(hit1)
                        for hit2 in region2.locations:
                            if hit1.query_id == hit2.query_id and hit1 != hit2:
                                hit2.identity = region2.identity
                                combine.append(hit2)
                        compare.append(combine)
                    
                    # Search through paired hits and fragments to find overlaps
                    # that can be removed
                    check_remove = []
                    addToRemoveList = True # Only remove a hit if this stays true
                    for pairs in compare:

                        #Dont remove unless all fragments have overlaps
                        if len(pairs) == 1:
                            addToRemoveList = False
                            break
                        elif len(pairs) == 2:

                            loc1 = pairs[0]
                            loc2 = pairs[1]

                            if loc1.reference_id == loc2.reference_id:
                                continue

                            overlap = min(loc1.query_stop, loc2.query_stop) - \
                            max(loc1.query_start, loc2.query_start) + 1
                        
                            overlap = max(0, overlap)

                            hit1_length = loc1.query_stop - loc1.query_start + 1
                            hit2_length = loc2.query_stop - loc2.query_start + 1
                            
                            if overlap >= min_overlap * min(hit1_length, hit2_length):
                                if loc1.identity > loc2.identity:
                                    check_remove.append(loc2.reference_id)
                                elif loc1.identity < loc2.identity:
                                    check_remove.append(loc1.reference_id)
                                else:
                                    if len(check_remove) == 0:
                                        check_remove.append(loc2.reference_id)
                        
                    # Check if all pairs agree that the same ref should be removed
                    # Then add to remove list
                    for ref1 in check_remove:
                        for ref2 in check_remove:
                            if ref1 != ref2:
                                addToRemoveList = False
                                break
                    
                    if addToRemoveList and len(check_remove) > 0:
                        to_remove.add(check_remove[0])

        for rm in to_remove:
            del best_hits[rm]

        for geno_region in best_hits.values():
            locus = geno_region.locations[0].reference_id.split("_")[0]
            allele = geno_region.locations[0].reference_id.split("_")[1].split("|")[0] ## Updated parsing for pathotypefinder
            gene_name = locus
            if geno_region.coverage >= blast_coverage:   
                results_out['results'][gene_name] = True ## True
                hit_information =[
                    {
                        'locus': locus,
                        'identity': geno_region.identity,
                        'allele': allele,
                        'hits': [
                            {
                                'contig_id': hit.query_id,
                                'locus':hit.reference_id,
                                'query_start': hit.query_start,
                                'query_stop': hit.query_stop,
                                'reference_start': hit.reference_start,
                                'reference_stop': hit.reference_stop,
                                'full_match': hit.full_match
                            } for hit in geno_region.locations
                        ]
                    }
                ]
                results_out['extra'].extend(hit_information)
            else:
                results_out['results'][gene_name] = False ## False

        for sequence_info in sequences.items():
            locus = sequence_info[0].split("_")[0]
            gene_name = locus
            if gene_name not in results_out['results']:
                results_out['results'][gene_name] = False ## False
  
        return results_out

    #pathotypefinder specific
    def get_pathotypes(self, loci):
        if any(locus.startswith('ipaH') for locus in loci) or 'ipaD' in loci:
            return ['EIEC/Shigella']

        if any(locus.startswith('stx') for locus in loci) and \
            ('aaiC' in loci or 'aggR' in loci or 'aatA' in loci or 'aap' in loci):
            return ['STEC/EAEC']

        if any(locus.startswith('stx') for locus in loci):
            return ['STEC']

        if 'ltcA' in loci or 'sta1' in loci or 'stb' in loci:
            return ['ETEC']

        if 'eae' in loci:
            if ('bfpA' in loci or 'eaf' in loci):
                return ['EPEC (typical)']
            else:
                return ['EPEC']

        if 'aaiC' in loci or 'aggR' in loci or 'aatA' in loci or 'aap' in loci:
            return ['EAEC']

        if 'daaC' in loci:
            return ['DAEC']

        return []

