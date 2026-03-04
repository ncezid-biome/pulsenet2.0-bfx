from collections import defaultdict
from csv import DictReader
from os import getenv
import os
from pathlib import Path
from shutil import move, rmtree
import json
from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.base.inputs import OrganismInput

from src.inputs import SeqSeroInputs
from src.outputs import SeqSeroOutputs
from src.config import CONTINGENCY_TABLES, LOOKUP_TABLES, MLST_VALIDATED


CONDA_BIN_DIR = os.getenv("CONDA_BIN_DIR", "/opt/conda/condabin") #Container
_QC_SECTION = "ANI"

_PULSENET_ORGANISM_ABBREVS = {
    "Salmonella": "SALM",
    "Escherichia": "STEC",
    "Listeria" : "LISTERIA",
    "Vibrio" : "VIBRIO",
    "Campylobacter" : "CAMPY",
    "Cronobacter" : "CRONO"
}


class SeqSero(Algorithm[SeqSeroInputs, SeqSeroOutputs]):

    outputs_class = SeqSeroOutputs

    def execute_stub(self):
        pass

    def execute_implementation(self):
        self.run_seqsero()
        seqsero_result = self.process_seqsero_out("SeqSero_result.txt")
        lookup_table = LOOKUP_TABLES[self.inputs.organism.genus]
        contingency_table = CONTINGENCY_TABLES[self.inputs.organism.genus]
        ispcr_file_contents = self.read_isPRC_out(self.inputs.ispcr_out)
        mlst_val_table = MLST_VALIDATED[self.inputs.organism.genus]
        mlst_st = self.read_MLST_out(self.inputs.mlst_out)

        results = {
        'formula': '',
        'serotype': ''
        }

        with open(str(self.inputs.species_ani_out), "r") as file: 
            ani_dict = json.load(file)
        
        if ani_dict["Subspecies"] == "":
            self.logger.info('No ANI value determined. Exiting salmonella serotyper.')
            results['serotype'] = ''
        else:
            ani_result = ani_dict["Subspecies"]
            results = self.determine_serotype(seqsero_result, ispcr_file_contents, contingency_table, ani_result, lookup_table, mlst_st, mlst_val_table)

        results_out = {
        'results': results,
        'extra': []
        }
        self.outputs.serotype.content = results_out


    def run_seqsero(self):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        cmd = ["mamba", "run"]
        cmd += ["-n", "seqsero"]
        cmd += ["SeqSero2s.py"]
        cmd += ["-m", "a"]
        cmd += ["-p", str(self.inputs.n_threads)]
        cmd += ["-t", "2"]
        cmd += ["-d", str(self.inputs.publish_dir)]
        cmd += ["-i", str(self.inputs.read1), str(self.inputs.read2)]
        process_stdout, err = run_external(cmd, self.logger, True, working_dir)
        self.logger.info(process_stdout)
        

    def process_seqsero_out(self, seqsero_file: Path) -> str:
        with open(seqsero_file) as fin:
            for line in fin:
                if line.startswith('Predicted antigenic profile:'):
                    return line.split('\t')[1].strip()

    def read_isPRC_out(self, pcr_file: Path) -> dict[str, bool]:
        self.logger.info(f"Reading isPCR output from {pcr_file}")
        with open(pcr_file) as fin:
            ispcr_output = json.load(fin)['results']
    
        return ispcr_output

    def read_MLST_out(self, mlst_file: Path) -> str:
        self.logger.info(f"Reading MLST ST output from {mlst_file}")
        with open(mlst_file) as fin:
            mlst_json = json.load(fin)

            try:
                mlst_st = mlst_json["result"]["mlst_results"][0]["sequence_type"]
            except:
                mlst_st = "-"
            
        return mlst_st

    def determine_serotype(
        self,
        seqsero_result: str,
        ispcr_file_contents: dict[str, bool],
        contingency_table: dict,
        ani_result: str,
        sslookuptable: dict,
        mlst_st: str,
        mlst_val_table: dict,
    ) -> str:
        results = {'formula': '', 'serotype': ''}
        formula = seqsero_result
        serotype = ''
        results['formula'] = '{spp} {formula}'.format(spp=ani_result,formula=formula)

        if sslookuptable.get(ani_result, False) and \
                sslookuptable[ani_result].get(seqsero_result, False):
            serotype = sslookuptable[ani_result][seqsero_result]
        else:
            serotype = '{spp} {formula}'.format(spp=ani_result,formula=formula)

        if serotype == "to genotyper":
            ispcr_contigency = self.check_ispcr_contingency(ispcr_file_contents, contingency_table, results['formula'])
            if ispcr_contigency == "Needs further review":
                results['serotype'] = '{spp} {formula}'.format(spp=ani_result,formula=formula)
            else:
                results['serotype'] = ispcr_contigency
        elif serotype == "to mlst":
            try:
                results['serotype'] = mlst_val_table[results['formula']][mlst_st]
            except:
                results['serotype'] = '{spp} {formula}'.format(spp=ani_result,formula=formula)
        else:
            results['serotype'] = serotype

        return results

    def check_ispcr_contingency(self, ispcr_results: dict[str, bool], run_contingency_table: dict, code: str) -> str:
        if code not in run_contingency_table:
            self.logger.info("serotype code not found in in silico PCR lookup table")
            return "Needs further review"
            
        for entry in run_contingency_table[code]:
            if all(contingency_result == ispcr_results[locus] for locus, contingency_result in entry.items() if locus != "BN Serotype"):
                self.logger.info(f"Found matching serotype based on in silico PCR result: {entry['BN Serotype']}")
                return entry["BN Serotype"]
        
        self.logger.info("Cound not find matching serotype based on in silico PCR result")
        return "Needs further review"