from collections import defaultdict
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

from src.inputs import EcKMASerotyperInputs
from src.outputs import EcKMASerotyperOutputs


CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")
_QC_SECTION = "STEC"


class EcKMASerotyper(Algorithm[EcKMASerotyperInputs, EcKMASerotyperOutputs]):

    outputs_class = EcKMASerotyperOutputs

    def execute_stub(self):
        working_dir = Path("tmp")
        working_dir.mkdir(exist_ok=True)
        cmd = [f"{CONDA_BIN_DIR}/kma"]
        process_stdout, _ = run_external(cmd, self.logger, True, working_dir)
        print(process_stdout)

    def execute_implementation(self):
        kdb_path = self.run_kma_index()
        self.kma_result_out = self.run_kma()
        self.parsed_result_out = self.parse_kma_output()
        ec_kma_serotyper_result_out = self.translate_results()
        self.outputs.ec_kma_serotyper_result_out.to_file()
        self.set_results()

    def run_kma_index(self):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        cmd = ["kma_index"]
        cmd += ["-i", str(self.inputs.otypes), str(self.inputs.htypes)]
        cmd += ["-o", "kma_db"]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return working_dir / "kma_db"

    def run_kma(self):
        #
        # https://bitbucket.org/genomicepidemiology/mykmafinder/src/master/
        #
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        kdb = self.run_kma_index()
        self.run_kma_index()
        cmd = ["kma"]
        cmd += ["-ipe", str(self.inputs.read1), str(self.inputs.read2)]
        cmd += ["-o", "ecoli_serotype"]
        cmd += ["-t_db", kdb]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return working_dir / "ecoli_serotype.res"

    def parse_kma_output(self):

        def check_antigen(template_hit):
            antigen = str(template_hit.split('_')[-1])
            if antigen.startswith('O'):
                return 'O'
            if antigen.startswith('H'):
                return 'H'

        def check_otype_groups(hits):
            # Check if valid O-type hits are in a related group. At this step,
            # we assume that each hit has passed the minimum parameter thresholds.
            pairings = [
                ['O62', 'O68'], ['O107', 'O117'], ['O90', 'O127'], ['O46', 'O134'],
                ['O20', 'O137'], ['O124', 'O164'], ['O153/O178', 'O153', 'O178'],
                ['O169/O183', 'O169', 'O183'], ['O123/O186', 'O123', 'O186'],
                ['O28/O42', 'O28', 'O42'], ['O50/O2', 'O2', 'O50'], ['O89', 'O101', 'O162'],
                ['O118/O151', 'O118', 'O151'], ['O13/O135', 'O13', 'O129', 'O135'],
                ['O17/O77', 'O17/O44', 'O17', 'O44', 'O73', 'O77', 'O106']
            ]
            checkpair = None
            for hit in hits:
                otype = hit['#Template'].strip().split('_')[-1]
                if checkpair is not None:
                    if otype not in checkpair:
                        # if passing O-type results are not in related groups
                        return False
                else:
                    for pair in pairings:
                        if otype in pair:
                            checkpair = pair
                            continue
            # if all passing O-type results are in related o-groups
            return True

        results_out = {
            'results': {
                'O': [],
                'H': []
            },
            'extra': []
        }

        with open(self.kma_result_out) as fh:
            records = DictReader(fh, delimiter="\t")
            headers = records.fieldnames
            if records:
                valid_hits = {"O": [], "H": []}
                for row in records:
                    template = row["#Template"].strip()
                    temp_identity = float(row["Template_Identity"].strip())
                    temp_coverage = float(row["Template_Coverage"].strip())
                    depth = float(row["Depth"].strip())

                    if (temp_identity >= float(60) and temp_coverage >= float(60) and depth >= float(5)):
                        hit_info = {}
                        for colname in headers:
                            hit_info[colname] = row[colname].strip()
                        valid_hits[check_antigen(template)].append(hit_info)

                # traverse the antigens once, filter for the top hit
                for antigen in valid_hits.keys():
                    hits = valid_hits[antigen]
                    if len(hits) == 0:
                        continue
                    elif len(hits) == 1:
                        top_hit = hits[0]
                        results_out["results"][check_antigen(top_hit["#Template"])].append(top_hit["#Template"].split("_")[-1])
                        results_out["extra"].append(top_hit)
                    else:
                        if antigen == "O":
                            if check_otype_groups(hits):
                                # if len(hits) > 1: returns the hit info (row) of the best hit per the Template Identity metric as a dict()
                                top_hit = max(hits, key=lambda x: x['Template_Identity'])
                                # use the parsed hit infor to populate the final result out dict
                                results_out['results'][check_antigen(top_hit['#Template'])].append(top_hit['#Template'].split('_')[-1])
                                results_out['extra'].append(top_hit)
                            else:
                                results_out['results']['O'].append('Needs Further Review')
                                # store the passing but non-related O-type hits for the SME to review
                                results_out['extra'].extend(hits)
                        else:
                            # handle H antigen results
                            # catch situations where 2 or more H-types score >= 99.0% template_id
                            # return dual H-types if conditions are met
                            count = 0
                            h_genes = []
                            h_hits = []
                            # first condition if >2 H-types with > 99% identity
                            for hit in hits:
                                if float(hit['Template_Identity']) >= float(99.0):
                                    count += 1
                                    gene_name = str(hit['#Template'].split('_')[0].strip())
                                    h_genes.append(gene_name)
                                    h_hits.append(hit)
                                else:
                                    continue
                            if count > 1:
                                dual_h_qualifiers = ['fllA', 'flkA', 'flmA']
                                # second condition fliC must be one of the top hits
                                if 'fliC' not in h_genes:
                                    # if fliC is not present but multiple high scoring H-types are still captured reports top hit
                                    top_hit = max(hits, key=lambda x: float(x['Template_Identity']))
                                    results_out['results'][check_antigen(top_hit['#Template'])].append(
                                        top_hit['#Template'].split('_')[-1])
                                    results_out['extra'].append(top_hit)
                                else:
                                    # final condition one of the qualifying rare H-genes must exist along side fliC
                                    check = any(h in dual_h_qualifiers for h in h_genes)
                                    if check is True:
                                        # report the combination of qualifying H types
                                        results_out['results']['H'].append(
                                            '/'.join([hit['#Template'].split('_')[-1] for hit in hits])
                                        )
                                    # add the details of the competing hits to the extra output for internal review
                                    for hit in h_hits:
                                        results_out['extra'].append(hit)
                            else:
                                # if there is a clear winner among valid hits
                                top_hit = max(hits, key=lambda x: float(x['Template_Identity']))
                                results_out['results'][check_antigen(top_hit['#Template'])].append(
                                    top_hit['#Template'].split('_')[-1])
                                results_out['extra'].append(top_hit)

                return results_out

            else:
                # return empty result dict if no results found
                return results_out

    def translate_results(self):
        # check for ambiguous O-types and replace with reportable O-types
        lookup_dict = {}
        with open(self.inputs.lookuptable, 'r') as fh:
            r = csv.reader(fh, delimiter="\t")
            next(r, None)
            for row in r:
                o_type = row[0].strip()
                o_report = row[1].strip()
                lookup_dict[o_type] = o_report

        # translate O-type to reportable format
        if self.parsed_result_out['results']['O']:
            o_type = self.parsed_result_out['results']['O'][0]
            if o_type in lookup_dict.keys():
                # replace with reportable o-type
                new_o_type = lookup_dict[o_type]
                self.parsed_result_out['results']['O'][0] = new_o_type
                self.outputs.ec_kma_serotyper_result_out.content = self.parsed_result_out
            else:
                self.outputs.ec_kma_serotyper_result_out.content = self.parsed_result_out
        else:
            # no O-types detected, skip translation
            self.outputs.ec_kma_serotyper_result_out.content = self.parsed_result_out

    def set_results(self):
        self.result = {
            "ec_kma_serotyper_result_out": self.inputs.publish_dir + str(self.outputs.ec_kma_serotyper_result_out.path)
        }
