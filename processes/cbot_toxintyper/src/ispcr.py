from collections import defaultdict, Counter
from csv import DictReader
from os import getenv
from os.path import exists
from pathlib import Path
from shutil import move, rmtree
from itertools import combinations
import json

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.base.inputs import OrganismInput
from ngs_pipeline_lib.tools.tools import gunzip_file
from src.inputs import ISPCRInputs
from src.outputs import ISPCROutputs
from src.ispcr_config import NUC_SIBLINGS, PCR_SETTINGS, PRIMERS ##
from src.genotype_classes import GenotypeHit, GenotypeResults


CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")

class ISPCR(Algorithm[ISPCRInputs, ISPCROutputs]):

    outputs_class = ISPCROutputs

    def execute_stub(self):
        pass

    def execute_implementation(self):
        self.working_dir = Path("./") #original
        organism = self.inputs.organism.genus
        self.pcr_settings = PCR_SETTINGS.get(organism)
        self.primers = PRIMERS.get(organism)
        
        self.assembly = gunzip_file(self.inputs.assembly)
        
        self.assembly_dict = {}
        with open(self.assembly) as reader:
            for header, seq in SimpleFastaParser(reader):
                self.assembly_dict[header] = seq

        # Add lookup by primer ID to primer dict
        temp_dict = {}
        for k,v in self.primers.items():
            for subdict in v:
                temp_dict[subdict["forward_primer_id"]] = v
                temp_dict[subdict["reverse_primer_id"]] = v
                temp_dict[subdict["id"]] = v

        self.primer_file = self.working_dir / "primers.fna"
        self.logger.info(f"Writing primers to {self.primer_file}")
        self.write_primers_fasta(self.primer_file)
        self.primers |= temp_dict

        self.blast_out = self.working_dir / "blast_out.txt"
        self.run_blast()
        self.process_blast()
        results = self.identify_amplicons()
        self.outputs.insilicopcr.content = results
        self.set_results()
        self.outputs.insilicopcr.to_file()


    def run_blast(self) -> None:
        #cmd = ["mamba", "run"]
        cmd = ["blastn"]
        cmd += ["-query", str(self.primer_file)]
        cmd += ["-subject", str(self.assembly)]
        cmd += ["-task", "blastn-short"]
        cmd += ["-out", str(self.blast_out)]
        cmd += ["-perc_identity", str(100.0*self.pcr_settings["PCR_PERCENT_IDENTITY"])]
        cmd += ["-outfmt", "6 std qlen"]
        cmd += ["-dust", "no"]
        blast_process_stdout, _ = run_external(
            cmd,
            self.logger,
            use_mamba_env=False,
            text=True,
            working_dir=self.working_dir
            )
        self.logger.info(blast_process_stdout)

        self.blast_results = GenotypeResults(str(self.blast_out))

        with open(str(self.blast_out), 'r+') as blastfile:
            contents = blastfile.read()
            blastfile.seek(0)
            blastfile.write(f"Query\tSubject\tPercent Identity\tAlignment Length\tMismatches\tGap Opens\tQuery Start\tQuery End\tSubject Start\tSubject End\tEvalue\tBit Score\tQuery Length\n{contents}")


    def process_blast(self) -> None:

        self.primer_hits = defaultdict(list)

        for hit in self.blast_results.hits:
            self.pad_hit(hit)
            
            # Adjust primer hit start and stop to be full length now the assembly hit is padded
            hit.query_start = 0
            hit.query_stop = hit.query_len - 1
            hit.query_seq = self.get_primer_seq(hit.query_id)

            hit.num_mismatches = self.check_mismatches(hit.query_seq, hit.reference_seq)

            hit.identity = 1.0 - (float(hit.num_mismatches) / float(hit.query_len))
            primer_id = self.primers[hit.query_id][0]["id"]
            self.primer_hits[primer_id].append(hit)

    def pad_hit(self, hit: GenotypeHit) -> None:
        contig = self.assembly_dict[hit.reference_id]
        if hit.forward:

            hit.reference_start = hit.reference_start - hit.query_start
            hit.reference_stop = hit.reference_stop + hit.query_len - hit.query_stop - 1

            pre_n = 'N' * (-1 * min(0, hit.reference_start))
            post_n = 'N' * (-1 * min(0, len(contig)-1 - hit.reference_stop))
            
            # It might be the case that the primer is hanging off the beginning
            # of a contig XOR the end.
            hit.reference_seq = pre_n + \
                contig[max(0, hit.reference_start):min(hit.reference_stop, len(contig)-1)+1] + \
                post_n
        
        else:

                hit.reference_start = hit.reference_start - \
                    hit.query_len + hit.query_stop + 1

                hit.reference_stop = hit.reference_stop + hit.query_start

                pre_n = 'N' * (-1 * min(0, len(contig)-1 - hit.reference_stop))
                post_n = 'N' * (-1 * min(0, hit.reference_start))

                new_r_seq = Seq(post_n + \
                    contig[max(0, hit.reference_start):min(hit.reference_stop, len(contig)-1)+1] + \
                    pre_n)
                hit.reference_seq = str(new_r_seq.reverse_complement())

    def get_primer_seq(self, primer: str) -> str:
        """ As primers dict is not set up to look up primer sequences directly,
        do it in a quick function to be tidier
        """
        # First get relevant dict from overall primers dict
        primer_sets = self.primers[primer]
        for primer_set in primer_sets:
            if primer == primer_set["forward_primer_id"]:
                return primer_set["forward_sequence"]
            elif primer == primer_set["reverse_primer_id"]:
                return primer_set["reverse_sequence"]
        
        return ''

    def check_mismatches(self, seq1: str, seq2: str) -> int:
        mismatches = 0

        for nuc1, nuc2 in zip(seq1, seq2):

            if nuc1 == nuc2:
                continue

            if nuc2 not in NUC_SIBLINGS[nuc1] \
                and nuc1 not in NUC_SIBLINGS[nuc2]:
                
                mismatches += 1

        return mismatches


    def identify_amplicons(self) -> None:

        # Store primer_hits that pass filtering
        good_primer_hits = {}

        self.logger.info(f"Before filtering, found {len(self.primer_hits)} potential primer amplicons")

        for primer_hit, hits in self.primer_hits.items():
            
            if len(hits) < 2:
                # need at least 2 primers to anneal (i.e., blast hits) for PCR
                continue

            # remove bad hits

            to_remove = []

            for i, hit in enumerate(hits):
                # count non-ATCG nucelotides and remove if over limit
                nuc_counts = Counter(hit.reference_seq)
                nonIUPAC = len(hit.reference_seq) - (
                    nuc_counts['A']
                    + nuc_counts['C']
                    + nuc_counts['T'] 
                    + nuc_counts['G']
                    )
                if nonIUPAC > self.pcr_settings["PCR_MAX_NONIUPAC"] or hit.num_mismatches > self.pcr_settings["PCR_MAX_MISMATCHES"]:
                    to_remove.append(i)

            for i in reversed(to_remove):
                del hits[i]
            
            if len(hits) >= 2:
                # If PCR amplicon still possible, keep the primer_hit
                good_primer_hits[primer_hit] = hits

        self.logger.info(f"After filtering, {len(good_primer_hits)} potential primer amplicons remain")

        # Next check if a reasonable primer pair exists for each set
        results = {'results': {subd["locus"]: False for v in self.primers.values() for subd in v}, "extra": []}
        for primer_hit, hits in good_primer_hits.items():
            for hit1, hit2 in combinations(hits, 2):
                
                # Skip if they don't hit the same contig
                if hit1.reference_id != hit2.reference_id:
                    continue

                # Skip if they are both the same direction (e.g., F and F primer)
                if hit1.forward == hit2.forward:
                    continue

                target_contig = self.assembly_dict[hit1.reference_id]

                # Get the length of the PCR product
                start = min(hit1.reference_start, hit2.reference_start)
                stop = max(hit1.reference_stop, hit2.reference_stop)

                length =  int(stop - start + 1)
                sequence = target_contig[start:stop+1]
                target_lengths = [primer_set["expected_length"] for primer_set in self.primers[primer_hit]]
                
                for target_length in target_lengths:
                    # If this amplicon differs from desired size by too much then skip
                    if not (1.0 - self.pcr_settings["PCR_MAX_LENGTH_DEVIATION"]) * target_length <= length <= \
                        (1.0 + self.pcr_settings["PCR_MAX_LENGTH_DEVIATION"]) * target_length:

                        continue
                    matched_target_length = target_length
                    if hit1.query_id in [subd["forward_primer_id"] for subd in self.primers[hit1.query_id]]:
                        forward = hit1
                        reverse = hit2
                    else:
                        forward = hit2
                        reverse = hit1

                    isforward = forward.reference_start < reverse.reference_start
                    primer_record = [x for x in self.primers[forward.query_id] if x["expected_length"] == matched_target_length][0]

                    result = {
                        'primer_id' : primer_hit,
                        'locus': primer_record["locus"],
                        'forward_primer_id': forward.query_id,
                        'forward_sequence': forward.query_seq,
                        'reverse_primer_id': reverse.query_id,
                        'reverse_sequence': reverse.query_seq,
                        'forward_mismatch': forward.num_mismatches,
                        'reverse_mismatch': reverse.num_mismatches,
                        'expected_len' : primer_record["expected_length"],
                        'actual_len': length,
                        'contig_id' : forward.reference_id,
                        'amplicon' : sequence,
                        'query_start': start,
                        'query_stop': stop,
                        'forward': isforward
                    }
                    results["results"][result["locus"]] = True
                    results["extra"].append(result)

        self.logger.info(f"in silico PCR results in {len([k for k,v in results['results'].items() if v])} amplicons of the expected length: "
                    + f"{', '.join([k for k,v in results['results'].items() if v])}")
        return results
    
    def write_primers_fasta(self, outfile: Path) -> Path:
        outcontents = ""
        written_primers = set()
        for id, ls in self.primers.items():
            for subd in ls:
                if (subd['forward_primer_id'], subd['forward_sequence']) not in written_primers:
                    outcontents += f">{subd['forward_primer_id']}\n{subd['forward_sequence']}\n"
                    written_primers.add((subd['forward_primer_id'], subd['forward_sequence']))
                if (subd['reverse_primer_id'], subd['reverse_sequence']) not in written_primers:
                    outcontents += f">{subd['reverse_primer_id']}\n{subd['reverse_sequence']}\n"
                    written_primers.add((subd['reverse_primer_id'], subd['reverse_sequence']))

        with open(outfile, 'w') as fout:
            fout.write(outcontents)
    

           
    def apply_quality_control(self, metrics, organism):
        quality_control = QualityControl(
            qc_dict=self.inputs.qc.get_dict(),
            report=self.qc_report,
            organism=organism
        )
        quality_control.report.add_metrics(metrics=metrics)
        quality_control.apply(section_name=_QC_SECTION, observations=metrics)

    def get_metrics(self):
        metrics = dict()
        metrics["Genus"] = self.outputs.best_hit.content["Genus"]
        metrics["Discrimination"] = self.outputs.best_hit.content["Discrimination"]
        metrics["Percent_Aligned"] = self.outputs.best_hit.content["Percent_Aligned"]
        metrics["ANI_Score"] = self.outputs.best_hit.content["ANI_Score"]

        return metrics

    def set_results(self):
        self.result = {
            "insilicopcr": self.inputs.publish_dir
            + str(self.outputs.insilicopcr.path)
        }
