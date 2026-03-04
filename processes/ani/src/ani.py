from collections import defaultdict
from csv import DictReader
from os import getenv
from pathlib import Path
from shutil import move, rmtree
from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.base.inputs import OrganismInput

from src.inputs import ANIInputs
from src.outputs import ANIOutputs


CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")
_ANI_SETTINGS = {
    "ALL": {
        "PERCENT_IDENTITY": 80.0,
        "MIN_COVERAGE": 70.0,
        "DISCRIMINATION": 2.0,
    }
}
_PULSENET_ORGANISM_ABBREVS = {
    "Salmonella": "SALM",
    "Escherichia": "STEC",
    "Listeria" : "LISTERIA",
    "Vibrio" : "VIBRIO",
    "Campylobacter" : "CAMPY",
    "Cronobacter" : "CRONO",
    "Clostridium" : "CBOT",
    "Yersinia" : "YERSINIA",
    "Grimontia": "Grimontia",
    "Photobacterium": "Photobacterium"
}


class ANI(Algorithm[ANIInputs, ANIOutputs]):

    outputs_class = ANIOutputs

    def execute_stub(self):
        working_dir = Path("tmp")
        working_dir.mkdir(exist_ok=True)
        cmd = [f"{CONDA_BIN_DIR}/mamba", "run"]
        cmd = ["perl", "./ani-m_1.pl"]
        cmd += ["--query", str(self.inputs.assembly)]
        cmd += ["--reference", str(self.inputs.reference)]
        cmd += ["--nThreads", str(self.inputs.n_threads)]
        cmd += ["--resultsdir", "."]
        process_stdout, _ = run_external(cmd, self.logger, True, working_dir)

    def execute_implementation(self):
        species_file = self.run_ani()
        results = self.parse_ani_output(species_file=species_file)
        self.interpret_ani(results)

        # The BMX wrapper requires that organism info is specified using
        # their OrganismInput class, so we have to put the ANI result in
        # an instance of that class to use it for QC rule lookup
        organism = OrganismInput(
            genus=_PULSENET_ORGANISM_ABBREVS.get(self.outputs.best_hit.content["Genus"]),
            species=self.outputs.best_hit.content["Species"]
            )

        qc_metrics = self.get_metrics()
        self.apply_quality_control(
            metrics=qc_metrics,
            organism=organism
            )
        
        print(_PULSENET_ORGANISM_ABBREVS.get(self.outputs.best_hit.content["Genus"]))
        self.set_results()
        self.outputs.best_hit.to_file()


    def run_ani(self):
        working_dir = Path(f"{self.inputs.publish_dir}")
        working_dir.mkdir(exist_ok=True)
        cmd = [f"{CONDA_BIN_DIR}/mamba", "run"]
        cmd = ["perl", "/app/ani-m_1.pl"]
        cmd += ["--query", str(self.inputs.assembly)]
        cmd += ["--reference", str(self.inputs.reference)]
        cmd += ["--nThreads", str(self.inputs.n_threads)]
        cmd += ["--resultsdir", str(working_dir)]
        process_stdout, _ = run_external(cmd, self.logger, False, working_dir)
        self.logger.info(process_stdout)
        move(working_dir / "results/raw/out.tsv", working_dir / "out.tsv")
        rmtree(working_dir / "results/")
        return working_dir / "out.tsv"
    
    def parse_ani_output(self, species_file):
        results = []
        with open(species_file) as flobj:
            table = iter(flobj)
            # Skip the header, to be thorough it is:
            # query reference   percent-aligned     ani     genus   species subspecies serotype
            next(table)

            for line in table:

                line = line.strip().split('\t')

                ani_line = {
                    "reference": line[1],
                    "percent_aligned": float(line[2]),
                    "ani_score": float(line[3]),
                    "taxonomy": tuple(line[4:])
                    }

                results.append(ani_line)
        return results
    
    def interpret_ani(self, ani_results):
        run_ani_settings = _ANI_SETTINGS["ALL"]
        perc_identity = run_ani_settings["PERCENT_IDENTITY"]
        coverage = run_ani_settings["MIN_COVERAGE"]
        discrimination = run_ani_settings["DISCRIMINATION"] # not used

        final_results = []

        for ani_line in ani_results:
            if ani_line["ani_score"] > perc_identity and \
                ani_line["percent_aligned"] > coverage:

                final_results.append(ani_line)
        
        best_hit_dict = {"Reference": "", "Percent_Aligned": 0, "ANI_Score": 0, "Discrimination": 0, "Genus" : "", "Species": "", "Subspecies": ""}

        if not len(final_results):
            self.outputs.best_hit.content = best_hit_dict
            return

        if len(final_results) == 1:
            self.ANI_best_hit_JSON(best_hit_dict, final_results[0], 100)
            return final_results[0]

        # Make sure that the first and the second hit are discriminatory
        final_results.sort(key = lambda x: -x["ani_score"])

        best_hit = final_results[0]
        second_best = final_results[1]

        hit_discrim = best_hit["ani_score"] - second_best["ani_score"]
        self.ANI_best_hit_JSON(best_hit_dict, best_hit, hit_discrim)

    

    def ANI_best_hit_JSON(self, best_hit_dict, input_dict, result_discrim):
        best_hit_dict["Reference"] = str(input_dict["reference"])
        best_hit_dict["Percent_Aligned"] = float(input_dict["percent_aligned"])
        best_hit_dict["ANI_Score"] = float(input_dict["ani_score"])
        best_hit_dict["Genus"] = str(input_dict["taxonomy"][0])
        best_hit_dict["Species"] = str(input_dict["taxonomy"][1])
        best_hit_dict["Discrimination"] = result_discrim
        try: 
            best_hit_dict["Subspecies"] = str(input_dict["taxonomy"][2])
        except: 
            best_hit_dict["Subspecies"] = ""
        self.outputs.best_hit.content = best_hit_dict
           
    def apply_quality_control(self, metrics, organism):
        quality_control = QualityControl(
            qc_dict=self.inputs.qc.get_dict(),
            report=self.qc_report,
            organism=organism
        )
        quality_control.report.add_metrics(metrics=metrics)
        quality_control.apply(section_name=self.inputs.qc_section, observations=metrics)

    def get_metrics(self):
        metrics = dict()
        metrics["Genus"] = self.outputs.best_hit.content["Genus"]
        metrics["Discrimination"] = self.outputs.best_hit.content["Discrimination"]
        metrics["Percent_Aligned"] = self.outputs.best_hit.content["Percent_Aligned"]
        metrics["ANI_Score"] = self.outputs.best_hit.content["ANI_Score"]

        return metrics

    def set_results(self):
        self.result = {
            "best_hit": self.inputs.publish_dir
            + str(self.outputs.best_hit.path)
        }
