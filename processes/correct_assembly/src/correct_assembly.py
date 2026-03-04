import os
from dataclasses import asdict
from logging import Logger
from pathlib import Path

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control import QCResult, QualityControl
from ngs_pipeline_lib.tools.quality_control.assembly import compute_metrics
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.tools.stub import generate_fake_fasta_file
from ngs_pipeline_lib.tools.tools import gunzip_file

from src.assembly_cleaning import (
    CleanedAssemblyInfo,
    assembly_cleaning,
    get_organism_cleaning_settings,
)
from src.inputs import CorrectAssemblyInputs
from src.outputs import CorrectAssemblyOutputs

##CONDA_BIN_DIR = os.getenv("CONDA_BIN_DIR", "/opt/conda/condabin") #Container


class CorrectAssembly(Algorithm[CorrectAssemblyInputs, CorrectAssemblyOutputs]):
    outputs_class = CorrectAssemblyOutputs

    def set_results(self):
        self.result = {
            "corrected_assembly": self.inputs.publish_dir
            + str(self.outputs.corrected_assembly.path),
        }

    def execute_stub(self):
        self.outputs.corrected_assembly.content = generate_fake_fasta_file()
        self.outputs.corrected_assembly.to_file()

        assembly_infos = compute_metrics(self.outputs.corrected_assembly.path)
        self.qc_report.add_metrics(asdict(assembly_infos))

        self.set_results()

    def execute_implementation(self):
        assembly = gunzip_file(self.inputs.assembly)

        quality_control = QualityControl(
            qc_dict=self.inputs.qc_kb.qc.get_dict(),
            organism=self.inputs.organism,
            report=self.qc_report,
        )

        alignment = self.create_alignment(
            self.inputs.bowtie,
            self.inputs.bowtie_n_hits,
            assembly,
            self.inputs.read1,
            self.inputs.read2,
            self.inputs.n_threads,
            self.logger,
        )

        binarized_alignment = self.binarize_alignment(alignment, self.logger)
        self.index_alignment(binarized_alignment, self.logger)

        cleaning_settings = get_organism_cleaning_settings(
            self.inputs.cleaning_kb.settings.get_dict(),
            self.inputs.organism.genus,
            self.inputs.organism.species,
        )
        cleaned_assembly = assembly.with_stem(f"cleaned_{assembly.stem}")
        cleaned_alignment = alignment.with_stem(f"cleaned_{binarized_alignment.stem}")
        cleaned_assembly_info = assembly_cleaning(
            assembly,
            cleaned_assembly,
            binarized_alignment,
            cleaned_alignment,
            self.outputs.depth_contigs.path,
            cleaning_settings,
            self.inputs.n_threads,
            self.logger,
        )

        self.index_fasta(cleaned_assembly, self.logger)

        self.cleanup_files(assembly)

        self.index_alignment(cleaned_alignment, self.logger)

        self.correct_assembly(
            cleaned_assembly,
            cleaned_alignment,
            self.outputs.corrected_assembly.path,
            self.inputs.n_threads,
            self.logger,
        )

        self.cleanup_files(cleaned_assembly)

        ##assembly_infos = compute_metrics(self.outputs.corrected_assembly.path) #Old metrics (incorrect)
        ##self.qc_report.add_metrics(asdict(assembly_infos))
        self.qc_report.add_metrics(self.get_additonal_metrics(cleaned_assembly_info))
        quality_control.apply(
            section_name="correctedAssembly",
            observations=self.qc_report.metrics,
        )

        corrected_alignment = self.create_alignment(
            self.inputs.bowtie,
            self.inputs.bowtie_n_hits,
            self.outputs.corrected_assembly.path,
            self.inputs.read1,
            self.inputs.read2,
            self.inputs.n_threads,
            self.logger,
        )

        self.index_fasta(self.outputs.corrected_assembly.path, self.logger)

        self.sort_and_compress_alignment(
            corrected_alignment,
            self.outputs.corrected_assembly.path,
            self.outputs.corrected_alignment.path,
            self.logger,
        )
        self.cleanup_files(self.outputs.corrected_assembly.path)

        self.set_results()

    @staticmethod
    def create_alignment(
        use_bowtie: bool,
        bowtie_n_hits: None | int,
        assembly: Path,
        read1: Path,
        read2: Path,
        threads: int,
        logger: Logger,
    ) -> Path:
        """
        Creates a .sam alignment file using the assembly and the reads
        Uses either bowtie or bwa
        """
        alignment = assembly.with_suffix(".sam")
        if use_bowtie:
            run_external(
                [
                    ##f"{CONDA_BIN_DIR}/mamba", #Container
                    "mamba", #Scicomp
                    "run",
                    "--no-capture-output",
                    "bowtie2-build",
                    str(assembly),
                    str(assembly),
                ],
                logger,
            )
            command_line = [
                ##f"{CONDA_BIN_DIR}/mamba", #Container
                "mamba", #Scicomp
                "run",
                "--no-capture-output",
                "bowtie2",
            ]
            command_line += ["--local"]
            if bowtie_n_hits and bowtie_n_hits > 1:
                command_line += ["-k", str(bowtie_n_hits)]
            else:
                command_line += ["-a"]
            command_line += ["-X", "1500"]
            command_line += ["-x", str(assembly)]
            command_line += ["-1", str(read1)]
            command_line += ["-2", str(read2)]
            command_line += ["-S", str(alignment)]
            command_line += ["-p", str(threads)]
            bowtie_stdout, _ = run_external(command_line, logger)
            logger.info(bowtie_stdout)
        else:
            run_external(
                [
                    ##f"{CONDA_BIN_DIR}/mamba", #Container
                    "mamba", #Scicomp
                    "run",
                    "--no-capture-output",
                    "bwa",
                    "index",
                    str(assembly),
                ],
                logger,
            )
            command_line = [
                ##f"{CONDA_BIN_DIR}/mamba", #Container
                "mamba", #Scicomp
                "--no-banner", #Scicomp
                "run",
                "--no-capture-output",
                "bwa",
                "mem",
            ]
            command_line += ["-t", str(threads)]
            command_line += [str(assembly)]
            command_line += [str(read1)]
            command_line += [str(read2)]
            mapping_output, _ = run_external(command_line, logger)
            with open(alignment, "w", encoding="utf-8") as open_file:
                open_file.write(mapping_output)
        return alignment

    @staticmethod
    def binarize_alignment(mapping: Path, logger: Logger) -> Path:
        """
        .sam -> .bam
        """
        binarized_mapping = (mapping.parent / mapping.stem).with_suffix(".bam")
        binarization_output, _ = run_external(
            [
                ##f"{CONDA_BIN_DIR}/mamba", #Container
                "mamba", #Scicomp
                "--no-banner", #Scicomp
                "run",
                "--no-capture-output",
                "samtools",
                "sort",
                str(mapping),
            ],
            logger,
            text=False,
        )
        with open(binarized_mapping, "wb") as open_file:
            open_file.write(binarization_output)
        # the .sam is not needed anymore, let's remove it
        mapping.unlink()
        return binarized_mapping

    @staticmethod
    def index_fasta(assembly: Path, logger: Logger) -> None:
        """
        .fasta -> .fasta.fai
        """
        run_external(
            [
                ##f"{CONDA_BIN_DIR}/mamba", #Container
                "mamba", #Scicomp
                "run",
                "--no-capture-output",
                "samtools",
                "faidx",
                str(assembly),
            ],
            logger,
        )

    @staticmethod
    def index_alignment(alignment: Path, logger: Logger) -> None:
        """
        .cram -> .cram.crai
        .bam -> .bam.bai
        .sam.gz -> .sam.gz.bai
        """
        run_external(
            [
                ##f"{CONDA_BIN_DIR}/mamba", #Container
                "mamba", #Scicomp
                "run",
                "--no-capture-output",
                "samtools",
                "index",
                str(alignment),
            ],
            logger,
        )

    @staticmethod
    def correct_assembly(
        assembly: Path,
        alignment: Path,
        output: Path,
        n_threads: int,
        logger: Logger,
    ) -> Path:
        ##cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "--no-capture-output", "pilon"] #Container
        cmd = ["mamba", "run", "--no-capture-output", "pilon"] #Scicomp
        cmd += ["-Xms512m"]
        cmd += ["-XX:MaxRAMPercentage=80"]
        cmd += ["--genome", str(assembly)]
        cmd += ["--frags", str(alignment)]
        cmd += ["--threads", str(n_threads)]
        cmd += ["--changes"]
        cmd += ["--output", output.stem]
        pilon_process_stdout, _ = run_external(cmd, logger)
        logger.info(pilon_process_stdout)

    @staticmethod
    def sort_and_compress_alignment(
        alignment: Path, assembly: Path, output: Path, logger: Logger
    ) -> None:
        """
        .bam -> .cram
        """
        sorted_alignment = Path("sorted.bam")
        run_external(
            [
                ##f"{CONDA_BIN_DIR}/mamba", #Container
                "mamba", #Scicomp
                "run",
                "--no-capture-output",
                "samtools",
                "sort",
                "-o",
                str(sorted_alignment),
                str(alignment),
            ],
            logger,
        )
        cmd = [
            ##f"{CONDA_BIN_DIR}/mamba", #Container
            "mamba", #Scicomp
            "run",
            "--no-capture-output",
            "samtools",
            "view",
        ]
        cmd += ["-h"]
        cmd += ["-o", str(output)]
        cmd += ["--output-fmt", "CRAM"]
        cmd += ["--reference", str(assembly)]
        cmd += [str(sorted_alignment)]
        run_external(
            cmd,
            logger,
        )
        sorted_alignment.unlink()

    @staticmethod
    def cleanup_files(filepath: Path) -> None:
        """
        Deletes files that follow the pattern "filename.*"
        """
        generated_files = Path(filepath.parent).glob(filepath.name + ".*")
        for generated_file in generated_files:
            generated_file.unlink()

    @staticmethod
    def get_additonal_metrics(
        cleaned_assembly_info: CleanedAssemblyInfo,
    ) -> dict[str, dict[str, float | int]]:
        return {
            "initial_length": cleaned_assembly_info.initial_assembly_length,
            "initial_average_depth": cleaned_assembly_info.initial_average_depth,
            "initial_n_contigs": cleaned_assembly_info.initial_contigs_number,
            "initial_gc": cleaned_assembly_info.initial_gc_content,
            "initial_general_horizontal_coverage": cleaned_assembly_info.initial_general_horizontal_coverage,
            "initial_n50": cleaned_assembly_info.initial_n50,
            "cleaned_length": cleaned_assembly_info.cleaned_assembly_length,
            "cleaned_average_depth": cleaned_assembly_info.cleaned_average_depth,
            "cleaned_n_contigs": cleaned_assembly_info.cleaned_contigs_number,
            "cleaned_gc": cleaned_assembly_info.cleaned_gc_content,
            "cleaned_general_horizontal_coverage": cleaned_assembly_info.cleaned_general_horizontal_coverage,
            "cleaned_n50": cleaned_assembly_info.cleaned_n50,
        }
