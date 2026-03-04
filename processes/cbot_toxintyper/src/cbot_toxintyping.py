from collections import defaultdict, Counter, namedtuple
from csv import DictReader
from os import getenv, makedirs, listdir, rename
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
import re
import json

from Bio.Seq import Seq

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.base.inputs import OrganismInput
from ngs_pipeline_lib.tools.tools import gunzip_file, gzip_file

from src.inputs import CbotToxinTypingInputs, ISPCRInputs
from src.outputs import CbotToxinTypingOutputs
from src.config import SETTINGS
from src.ispcr_config import PCR_SETTINGS, PRIMERS
from src.ispcr import ISPCR
from src.genotype_classes import GenotypeHit, GenotypeResults, SequenceInfo


class CbotToxinTyping(Algorithm[CbotToxinTypingInputs, CbotToxinTypingOutputs]):

    outputs_class = CbotToxinTypingOutputs

    def execute_stub(self):
        pass

    def execute_implementation(self):
        self.working_dir = Path("./")
        blast_out = Path(f"{self.working_dir}_{self.inputs.sample_id}_cbottoxintyping_blast.txt")
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

        ## Get top hits from results
        top_hits = self.interpret_blast( 
            results,
            sequences,
            percent_identity = settings["PERCENT_IDENTITY"], 
            max_gaps = settings["MAX_GAPS"],
            min_alignment_percentage = settings["MIN_ALIGNMENT_PERCENTAGE"]
        )
        self.logger.info("Finished interpreting BLAST results")

        if not top_hits:
            self.logger.warn('No BLAST results were retained')
            empty_results = self.write_empty_results()

            self.outputs.cbot_toxintyper.content = empty_results
            self.outputs.cbot_toxintyper.to_file()

        else:
            #save top hits to individual fasta files to run through ispcr
            self.logger.info('Writing out results from initial BLAST...')
            self.write_individual_blasts(top_hits = top_hits)
            self.logger.info('Initial BlAST results written out')

            self.logger.info('Running in silico PCR...')

            #Grab PCR settings and primers
            self.pcr_settings = PCR_SETTINGS[genus]
            self.pcr_primers = PRIMERS[genus]

            #Run isPCR on the top blast hits
            pcr_results = self.run_insilico_pcr(
                top_hits = top_hits, 
            )
            self.logger.info('Finished running in silico PCR')
    
            #Write final results
            self.logger.info('Writing out top two results...')
            final_results = self.write_final_results(pcr_results)
            self.logger.info('Finished writing out top two results...')

            self.outputs.cbot_toxintyper.content = final_results
            self.outputs.cbot_toxintyper.to_file()


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

        cmd = ["blastn"]
        cmd += ["-query", str(query)]
        cmd += ["-db", db_path]
        cmd += ["-task", "blastn"]
        cmd += ["-out", str(outfile)]
        cmd += ["-perc_identity", str(100.0*percent_identity)]
        cmd += ["-outfmt", "6 std qlen qseq sseq"]
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

    #Interpret Cbot toxin blast results
    def interpret_blast(self, results:GenotypeResults, sequences:dict, percent_identity:float, max_gaps:float, min_alignment_percentage:float) -> list:

        self.logger.info("Interpreting BLAST results...")

        # Create a structure to help hold hit information and errors for scoring system
        Hit = namedtuple('Hit', ['contig', 'subtype', 'blast_hit', 'errors', 'gapless' ])
        top_hits = []
        
        for hit in results.hits:
            errors = {'0': False, '1': False, '2': False, '3': False, '4': False, '4.1': False, '4.2': False,
                 '4.3': False, '5': False, '5.1': False, '5.2': False , '5.3': False , '6': False}
            
            # Pass through only hits that meet our criteria
            if hit.num_gap_opens <= max_gaps and hit.identity >= percent_identity and hit.relative_len >= min_alignment_percentage:
                
                # Scoring 2^2
                #Check if blast hit contains gaps
                if hit.num_gap_opens > 0:
                    errors['2'] = True
                
                # Ungap the hit sequence and use biopython to reverse_complement it if the hit is not forward
                seq = Seq(hit.query_seq)
                gapless = seq.replace("-", "")

                if not hit.forward:
                    gapless = gapless.reverse_complement()
                
                # Use biopython to translate the hit sequence into protein for scoring system
                translated = gapless.translate()

                # Scoring 2^4
                #Check if translated blast hit does not start with M
                if translated[0] != 'M':
                    errors['4'] = True
                    errors['4.1'] = True

                # Check if translated blast hit has more then 2 stop codons
                if len(translated.split('*')) > 2:
                    errors['4'] = True
                    errors['4.2'] = True

                # Check if translated blast hit gaps are out of frame
                gap_positions = list(re.finditer('-+', str(seq)))                
                for gap in gap_positions:
                    if (gap.end() - gap.start()) % 3 != 0:
                        errors['4'] = True
                        errors['4.3'] = True
                        break
                    
                # Get a clean contig name
                contig_hit = hit.query_id.replace('|', '_') + '_'
                contig_hit = contig_hit.replace(',', '_')
                contig_hit = contig_hit.replace('/', '_')
    
                # Create and append our top hit structure with hit and error information
                top_hit = Hit(blast_hit=hit, contig=contig_hit, subtype=hit.reference_id, errors=errors, gapless=gapless)
                top_hits.append(top_hit)

        return top_hits

    def write_empty_results(self):
        
        # Write out results for no BLAST hits
        B_err = {
            "Missing_Blast_Hits": True,
            "Low_Percent_Hit_Length": False,
            "Multi_Region_Hit": False,
            "Multi_Contig_Hit": False,
            "Hit_Contains_Gaps": False,
            "Translated_Missing_Start": False,
            "Translated_Internal_Stops": False,
            "Translated_Gap_Out_Of_Frame": False,
            "isPCR_Results_Missing": False,
            "isPCR_BLAST_Discrepancy": False,
            "isPCR_Length_Discrepancy": False
        }
        BoNT_1 = {"BoNT_Seq": "", "Subtype": "NT", "Subtype_QC": "FAIL", "PercentHitLength": "", "BoNT_errors": B_err}
        
        results_out = {
            'results': {
                'ToxinSubtype': 'NT',
                'BoNT_1': BoNT_1,
            },
            'extra': []
        }     
        return results_out

    def write_individual_blasts(self, top_hits:list):
        
        for hit in top_hits:
            individual_path = join(self.working_dir, hit.contig + hit.subtype)

            if not exists(individual_path):
                makedirs(individual_path)
        
            file_path = join(individual_path, 'hit.fasta')
            
            with open(file_path, 'w') as fh:
                fh.write('>' + hit.contig + hit.subtype + '\n')
                fh.write(str(hit.gapless))

            gzip_file(Path(file_path))

    def run_insilico_pcr(self, top_hits:list):
        
        # Iterate through all the previously created fasta files and perform
        # insilico pcr
        for hit in top_hits:
            
            individual_path = join(self.working_dir, hit.contig + hit.subtype)
            file_path = join(individual_path, 'hit.fasta.gz')

            self.logger.info(f"Running insilico PCR for {hit.contig}")
            
            pcr_inputs = ISPCRInputs(
                sample_id = self.inputs.sample_id,
                organism = self.inputs.organism,
                logging_dir = f"{individual_path}/logs",
                publish_dir = individual_path,
                assembly = file_path,
            )

            pcr_data = ISPCR(pcr_inputs)
            pcr_data.execute_implementation()
            
            self.logger.info(f"Finished running insilico PCR for {hit.contig}")

            insilico_pcr_path = join(self.working_dir, 'insilicopcr.json')

            # Write out insilicopcr results for each top hit to interpret later
            ispcr_results_flag = 0 # if this stays zero then no pcr results found
            if exists(insilico_pcr_path):

                new_insilico_path = join(self.working_dir, f"{hit.contig + hit.subtype}.insilicopcr.json")
                rename(insilico_pcr_path, new_insilico_path)
            
                with open(new_insilico_path, 'r') as fh:
                    json_data = json.load(fh)

                if json_data is not None:
                    #Scoring 2^5
                    subtype = hit.subtype[0]
                    
                    subtype_list = []
                    for k, v in json_data["results"].items():
                        if v:
                            subtype_list.append(k)

                    #If not empty, then check expected vs actual lengths
                    if json_data["extra"]:
                        ispcr_results_flag = 1 #results has a true
                        expected_length = json_data["extra"][0]["expected_len"]
                        actual_length = json_data["extra"][0]["actual_len"]

                        # Check whether there are more than one PCR called, the called PCR is not the subtype,
                        # or whether the actual length is not in range of the expected_length
                        # Actual length deviates by 1bp in insilico PCR
                        if len(subtype_list) > 1 or subtype not in subtype_list:
                            self.logger.warn(f"PCR subtype hit and BLAST subtype hit have discrepancies for {hit.contig}")
                            hit.errors['5'] = True
                            hit.errors['5.1'] = True
                        elif actual_length not in [expected_length, expected_length + 1]:
                            self.logger.warn(f"PCR subtype hit length not in expected range for {hit.contig}")
                            hit.errors['5'] = True
                            hit.errors['5.2'] = True
            
            # Add this error if file is missing or pcr results are null or empty
            if ispcr_results_flag == 0:
                self.logger.warn(f"No insilico_pcr results were found for {hit.contig}")
                hit.errors['5'] = True
                hit.errors['5.3'] = True
   
        return top_hits

    def write_final_results(self, hits):
        
        results_out = {
            'results': {
                'ToxinSubtype': '',
            },
            'extra': []
        }

        # Sort hits by percent_id and alignment
        sorted_hits = sorted(hits, key=lambda hit: (hit.blast_hit.identity, hit.blast_hit.relative_len), reverse=True)

        toxinsubtype = {}
        prev_subtype = ''
        error_dict = {}

        # Iterate through sorted hits to keep top two
        for hit in sorted_hits:
            top_subtype = hit.subtype[0]       
            contig = hit.contig.split('_')[0]
            start = hit.blast_hit.query_start
            stop = hit.blast_hit.query_stop

            # Last minute error check 2^0 and 2^1:
            if top_subtype not in error_dict or error_dict[top_subtype]['0'] == 0 or error_dict[top_subtype]['1'] == 0:
            
                error_dict[top_subtype] = {'0': False, '1': False, 'contig_list': [contig], 'region_list': [(start, stop)]}
            
            else:
                #error 2^1
                if contig not in error_dict[top_subtype]['contig_list']:
                    error_dict[top_subtype]['1'] = True
                    error_dict[top_subtype]['contig_list'].append(contig)
                
                #error 2^0
                for region in error_dict[top_subtype]['region_list']:
                    if start >= region[1] or stop <= region[0]:
                        error_dict[top_subtype]['0'] = True
                    
                else:
                    error_dict[top_subtype]['region_list'].append((start, stop))      

            
            # Check whether we've seen the subtype letter before or if there are already two top hits
            if top_subtype not in toxinsubtype and top_subtype != prev_subtype and len(toxinsubtype) < 2:
                
                toxinsubtype[hit.subtype.split("_")[0]] = hit
                
                # Check to see if alignment is less than 100%
                if min([hit.blast_hit.relative_len, 1.0]) != 1.0:
                    hit.errors['6'] = True ##
            
                prev_subtype = top_subtype
        
        
        # Sort dictionary first by keys to print subtype results alphabetically
        toxin_keys = list(toxinsubtype.keys())
        toxin_keys.sort()
        sorted_toxinsubtype = {y: toxinsubtype[y] for y in toxin_keys}

        tox_silence = ""
        toxin_results = []
        i = 1
        for k, v in sorted_toxinsubtype.items():

            #Add new key
            results_out["results"][f"BoNT_{i}"] = {}

            # Save BLAST sequence of top hits
            results_out["results"][f"BoNT_{i}"]["BoNT_Seq"] = v.blast_hit.query_seq
            
            # Update the errors for the top hits
            v.errors['0'] = error_dict[k[0]]['0']
            v.errors['1'] = error_dict[k[0]]['1']

            #Print subtype
            results_out["results"][f"BoNT_{i}"][f"Subtype"] = k
            
            # Write warning level for QC
            results_out["results"][f"BoNT_{i}"][f"Subtype_QC"] = 'PASS'
            # Check for FAIL first
            if v.errors['1']:
                results_out["results"][f"BoNT_{i}"][f"Subtype_QC"] = 'FAIL'
            # If no FAIL, check for WARN
            elif v.errors['0'] or v.errors['2'] or v.errors['4'] or v.errors['5'] or v.errors['6']:
                results_out["results"][f"BoNT_{i}"][f"Subtype_QC"] = 'WARN'            

            # Write out percent aligned
            results_out["results"][f"BoNT_{i}"][f"PercentHitLength"] = round(min(100*v.blast_hit.relative_len, 100), 2)
            
            # Tag silent subtypes and non passing QC and 
            # add to subtype list
            if v.errors['4']:
                tox_silence = str(k)+"(silent)"
                toxin_results.append(tox_silence)
            #This will cover QC FAIL and WARN
            elif True in v.errors.values():
                tox_NFR = str(k)+"(Needs Further Review)"
                toxin_results.append(tox_NFR)
            else:
                toxin_results.append(k)


            ## Set error stats
            # 0 = subtype BLASTs to multiple regions on same contig
            # 1 = subtype BLASTs to multiple contigs
            # 2 = blast hit contains gaps
            # 3 = nothing - used to be fragments
            # 4 Translated blast hit problems
            # 4.1 = translated blast hit does not start with M, 
            # 4.2 = translated blast hit has more then 2 stop codons, 
            # 4.3 = translated blast hit gaps are out of frame
            # 5 isPCR problems
            # 5.1 = PCR blast discrepency
            # 5.2 = PCR length discrepency 
            # 5.3 = No insilico_pcr results were found
            # 6 = Low percent hit length set to true if < 100%
            
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"] = {}
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Missing_Blast_Hits"] = False
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Low_Percent_Hit_Length"] = v.errors['6']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Multi_Region_Hit"] = v.errors['0']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Multi_Contig_Hit"] = v.errors['1']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Hit_Contains_Gaps"] = v.errors['2']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Translated_Missing_Start"] = v.errors['4.1']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Translated_Internal_Stops"] = v.errors['4.2']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"Translated_Gap_Out_Of_Frame"] = v.errors['4.3']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"isPCR_Results_Missing"] = v.errors['5.3']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"isPCR_BLAST_Discrepancy"] = v.errors['5.1']
            results_out["results"][f"BoNT_{i}"]["BoNT_errors"][f"isPCR_Length_Discrepancy"] = v.errors['5.2']

            i += 1

        results_out["results"]["ToxinSubtype"] = ",".join(toxin_results)
        
        return results_out
