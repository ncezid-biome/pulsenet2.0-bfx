import os
from dataclasses import asdict
from pathlib import Path
from shutil import rmtree

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control import QualityControl
from ngs_pipeline_lib.tools.quality_control.assembly import compute_metrics
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.tools.stub import generate_fake_fasta_file

from src.inputs import GenerateAssemblyInputs
from src.outputs import GenerateAssemblyOutputs

CONDA_BIN_DIR = os.getenv("CONDA_BIN_DIR", "/opt/conda/condabin") #Container
SPADES_FOLDER = "spades"


class GenerateAssembly(Algorithm[GenerateAssemblyInputs, GenerateAssemblyOutputs]):
    outputs_class = GenerateAssemblyOutputs

    def execute_stub(self):
        self.outputs.assembly.content = generate_fake_fasta_file()
        self.outputs.assembly.to_file()

        assembly_infos = compute_metrics(self.outputs.assembly.path)
        self.qc_report.add_metrics(asdict(assembly_infos))

        self.set_results()

    def execute_implementation(self):
        self.call_spades()

        self.clean()

        quality_control = QualityControl(
            qc_dict=self.inputs.qc_kb.qc.get_dict(),
            organism=self.inputs.organism,
            report=self.qc_report,
        )
        assembly_infos = compute_metrics(self.outputs.assembly.path)
        assembly_infos.gc = round(assembly_infos.gc*100, 2) #fix for GC format
        self.qc_report.add_metrics(asdict(assembly_infos))

        quality_control.apply(
            section_name="rawAssembly", observations=self.qc_report.metrics
        )

        self.set_results()

    def set_results(self):
        self.result = {
            "assembly": self.inputs.publish_dir + str(self.outputs.assembly.path),
        }

    def call_spades(self):
        spades_folder = Path(SPADES_FOLDER)
        cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "spades.py"] #Container
        #cmd = ["spades.py"] # scicomp
        cmd += ["-1", str(self.inputs.read1)]
        cmd += ["-2", str(self.inputs.read2)]
        cmd += ["-o", str(spades_folder)]
        cmd += ["-t", str(self.inputs.n_threads)]
        cmd += ["--isolate"]
        run_external(cmd, self.logger)
        (spades_folder / "contigs.fasta").rename(self.outputs.assembly.path)

    def clean(self):
        spades_folder = Path(SPADES_FOLDER)
        rmtree(spades_folder)
