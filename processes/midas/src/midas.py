from collections import defaultdict
from csv import DictReader
from os import getenv
from pathlib import Path
from shutil import move
from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from subprocess import run, PIPE

from src.inputs import MidasInputs
from src.outputs import MidasOutputs


CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")
_QC_SECTION = "MIDAS"


class Midas(Algorithm[MidasInputs, MidasOutputs]):

    outputs_class = MidasOutputs

    def execute_stub(self):
        working_dir = Path("tmp")
        working_dir.mkdir(exist_ok=True)
        cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "run_midas.py"]
        process_stdout, _ = run_external(cmd, self.logger, True, working_dir)
        print(process_stdout)

    def execute_implementation(self):
        species_file = self.run_midas()
        self.parse_midas_output(species_file=species_file)
        qc_metrics = self.get_metrics()
        self.apply_quality_control(metrics=qc_metrics)
        self.outputs.midas_result.to_file()
        self.set_results()

    def run_midas(self):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        cmd = ["python2.7", "/opt/conda/envs/midas/bin/run_midas.py"]
        cmd += ["species", str(working_dir.resolve())]
        cmd += ["-1", str(self.inputs.read1)]
        cmd += ["-2", str(self.inputs.read2)]
        cmd += ["-t", str(self.inputs.n_threads)]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        return working_dir / "species/species_profile.txt"

    def parse_midas_output(self, species_file: Path):
        _valids = defaultdict(float) # never really used
        _species = defaultdict(float) # never really used
        _all_hits = []

        with open(species_file) as fh:
            reader = DictReader(fh, delimiter='\t')
            for result in reader:
                result_data = {}
                result_data['species_id'] = result.get('species_id')
                result_data['count_reads'] = float(result.get('count_reads'))
                result_data['coverage'] = float(result.get('coverage'))
                result_data['relative_abundance'] = float(result.get('relative_abundance'))

                genus_species = '_'.join(result_data['species_id'].split('_')[:2])

                _valids[genus_species] += result_data['coverage']
                _species[genus_species] += result_data['relative_abundance']
                _all_hits.append(result_data)
            # Abstract resulting dict() by genus name, hits with >= 0.0% coverage, and count occurrences
        genus_results_out = {'results': [], 'extra': []}

        s = {}
        for entry in _all_hits:
            if float(entry['coverage']) <= float(0):
                continue
            else:
                genus = entry['species_id'].split('_')[0]
                coverage = float(entry['coverage'])
                relative_abundance=float(entry['relative_abundance'])
                if genus not in s:
                    s[genus] = {'genus': genus, 'coverage': coverage, 'relative_abundance': relative_abundance, 'count': int(1)}
                else:
                    new_cov = s[genus]['coverage'] + entry['coverage']
                    new_ab = s[genus]['relative_abundance'] + entry['relative_abundance']
                    new_count = s[genus]['count'] + 1
                    s[genus].update({'coverage': new_cov,'relative_abundance':new_ab,'count': new_count})

        # All hits with > 0.00% coverage added to 'Extra' results
        top_species_hits = sorted(filter(lambda x: x['coverage'] > float(0), _all_hits), key=lambda x: x['coverage'], reverse=True)
        # Condensed results
        all_filtered_hits = [(k, v) for k, v in s.items()]
        sorted_filtered_hits = sorted(all_filtered_hits, key=lambda x: x[1]['coverage'], reverse=True)
        for k, v in sorted_filtered_hits:
            d = {}
            d[k] = v
            genus_results_out['results'].append(d[k])


        genus_results_out['extra'].extend(top_species_hits)
        self.outputs.midas_result.content = genus_results_out
       
    def apply_quality_control(self, metrics):
        quality_control = QualityControl(
            qc_dict=self.inputs.qc_kb.qc.get_dict(),
            report=self.qc_report,
            organism=self.inputs.organism,
        )
        quality_control.report.add_metrics(metrics=metrics)
        qc_metrics = {
            "primaryCoverage": metrics["primaryCoverage"]["coverage"],
            "secondaryCoverage": metrics["secondaryCoverage"]["coverage"] 
        }
        if metrics["primaryCoverage"]["coverage"] > 0:
            qc_metrics.update({"primaryGenus": metrics["primaryCoverage"]["genus"]})
        else:
            qc_metrics.update({"primaryGenus": ""})
        
        if metrics["secondaryCoverage"]["coverage"] > 0:
            qc_metrics.update({"secondaryGenus": metrics["secondaryCoverage"]["genus"]})
        else:
            qc_metrics.update({"secondaryGenus": ""})

        quality_control.apply(section_name=_QC_SECTION, observations=qc_metrics)

    def get_metrics(self):
        metrics = dict()
        #print genus to stdout 
        if len(self.outputs.midas_result.content["results"]) >= 1:
            print(self.outputs.midas_result.content["results"][0]["genus"], end="")
            metrics["primaryCoverage"] = self.outputs.midas_result.content["results"][0]
        else: 
            print("Unidentified", end="")
            metrics["primaryCoverage"] = {"coverage": 0}

        if len(self.outputs.midas_result.content["results"]) >= 2:
            metrics["secondaryCoverage"] = self.outputs.midas_result.content["results"][1]
        else:
            metrics["secondaryCoverage"] = {"coverage": 0}
        return metrics
    
    def set_results(self):

        self.result = {
            "midas_result_out": self.inputs.publish_dir + str(self.outputs.midas_result.path),
            "species_profile": self.inputs.publish_dir + str(self.outputs.species_profile.path),
        }
