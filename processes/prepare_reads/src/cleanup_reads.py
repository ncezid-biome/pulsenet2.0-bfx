import os
from gzip import open as gzip_open
from json import load
from logging import Logger
from pathlib import Path

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.base.inputs import OrganismInput
from ngs_pipeline_lib.tools.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.tools.stub import generate_fake_fastq_file

from src import stub
from src.inputs import CleanupReadsInputs
from src.outputs import CleanupReadsOutputs

CONDA_BIN_DIR = os.getenv("CONDA_BIN_DIR", "/opt/conda/condabin") #Container

class CleanupReads(Algorithm[CleanupReadsInputs, CleanupReadsOutputs]):
    outputs_class = CleanupReadsOutputs

    def set_results(self):

        self.result = {
            "cleaned_read_1": self.inputs.publish_dir
            + str(self.outputs.cleaned_read_1.path),
            "cleaned_read_2": self.inputs.publish_dir
            + str(self.outputs.cleaned_read_2.path),
        }

    def execute_stub(self):
        self.outputs.cleaned_read_1.content = generate_fake_fastq_file()
        self.outputs.cleaned_read_2.content = generate_fake_fastq_file()
        self.outputs.fastp_report.content = stub.FAKE_FASTP_REPORT

        self.outputs.cleaned_read_1.to_file()
        self.outputs.cleaned_read_2.to_file()
        self.outputs.fastp_report.to_file()

        self.set_results()

    def execute_implementation(self):
        self.check_corruption(self.inputs.read1)
        self.check_corruption(self.inputs.read2)

        downsampled_read_1 = Path("downsampled_read_1.fastq.gz")
        downsampled_read_2 = Path("downsampled_read_2.fastq.gz")

        quality_control = QualityControl(
            qc_dict=self.inputs.qc_kb.qc.get_dict(),
            organism=self.inputs.organism,
            report=self.qc_report,
        )
        genome_size = self.get_organism_value(
            self.inputs.genome_kb.sizes.get_dict(), self.inputs.organism
        )
        target_depth = self.get_organism_value(
            self.inputs.genome_kb.target_depths.get_dict(), self.inputs.organism, 150
        )

        initial_estimated_vertical_depth = self.sample_fastq(
            self.inputs.read1,
            self.inputs.read2,
            downsampled_read_1,
            downsampled_read_2,
            genome_size,
            target_depth,
            self.logger,
        )
        self.call_fastp(
            self.outputs.fastp_report.path,
            downsampled_read_1,
            downsampled_read_2,
            self.outputs.cleaned_read_1.path,
            self.outputs.cleaned_read_2.path,
            self.logger,
        )
        postSubsamplingEstimatedVerticalDepth = (
            self.get_post_subsampling_estimated_vertical_depth(
                self.outputs.fastp_report.path, genome_size
            )
        )
        fastp_observations = self.get_fastp_observations(self.outputs.fastp_report.path)
        fastp_observations[
            "initialEstimatedVerticalDepth"
        ] = initial_estimated_vertical_depth
        fastp_observations[
            "postSubsamplingEstimatedVerticalDepth"
        ] = postSubsamplingEstimatedVerticalDepth

        quality_control.apply(
            section_name="fastpThresholds", observations=fastp_observations
        )

        quality_control.report.add_metrics(metrics=fastp_observations)

        self.set_results()

    @staticmethod
    def check_corruption(compressed_file: Path) -> None:
        """
        Attemps to decompress and read the file, will raise an error if it fails
        """
        with gzip_open(str(compressed_file), "rt", encoding="utf-8") as open_file:
            while True:
                try:
                    next(open_file)
                except StopIteration:
                    return

    @staticmethod
    def get_organism_value(
        organism_to_value: dict[str, int], organism: OrganismInput, default: int = 0
    ) -> int:
        return organism_to_value.get(
            f"{organism.genus} {organism.species}",
            organism_to_value.get(
                organism.genus, organism_to_value.get("DEFAULT", default)
            ),
        )

    @staticmethod
    def get_post_subsampling_estimated_vertical_depth(
        fastp_report: Path, genome_size: int
    ) -> float:
        with open(fastp_report, encoding="utf-8") as reader:
            fastp_report_json = load(reader)
        total_bases = fastp_report_json["summary"]["after_filtering"]["total_bases"]
        try:
            return round(total_bases / (genome_size * 1e6), 2)
        except ZeroDivisionError:
            return float("nan")

    @staticmethod
    # FIXME Original depth lower than target depth ->action?
    def sample_fastq(
        read_1: Path,
        read_2: Path,
        output_filename_1: Path,
        output_filename_2: Path,
        genome_size: int,
        target_depth: int,
        logger: Logger,
    ) -> float:
        def get_bases(reads: Path, logger: Logger) -> int:
            seqtk_output, _ = run_external(
                [f"{CONDA_BIN_DIR}/mamba", "run", "seqtk", "fqchk", str(reads)], logger #Container
                #[f"seqtk", "fqchk", str(reads)], logger #Scicomp
            )
            return int(seqtk_output.splitlines()[2].split()[1])

        def sample_and_compress(
            reads: Path, ratio: float, output: Path, logger: Logger
        ):
            seqtk_output, _ = run_external(
                [
                    f"{CONDA_BIN_DIR}/mamba", #Container
                    "run", #Container
                    "seqtk",
                    "sample",
                    "-s100",
                    str(reads),
                    str(ratio),
                ],
                logger,
                text=False,
            )
            with gzip_open(str(output), "w") as output_file:
                output_file.write(seqtk_output)

        bases_1 = get_bases(read_1, logger)
        logger.info("Bases P1:" + str(bases_1))

        bases_2 = get_bases(read_2, logger)
        logger.info("Bases P2:" + str(bases_2))
        try:
            estimated_coverage = (bases_1 + bases_2) / (genome_size * 1e6)
            logger.info("Estimated coverage: " + str(estimated_coverage))
            ratio = target_depth / estimated_coverage

            initial_estimated_vertical_depth = round(estimated_coverage, 2)
        except ZeroDivisionError:
            initial_estimated_vertical_depth = float("nan")
            # FIXME check if it works with symlink
            output_filename_1.symlink_to(read_1)
            output_filename_2.symlink_to(read_2)
            return initial_estimated_vertical_depth
        else:
            logger.info("Subsample target ratio:" + str(ratio))
            if ratio < 1:
                logger.info(f"Writing {output_filename_1}")
                sample_and_compress(read_1, ratio, output_filename_1, logger)
                logger.info(f"Writing {output_filename_2}")
                sample_and_compress(read_2, ratio, output_filename_2, logger)
                logger.info("All done. Have a nice day!")
            else:
                output_filename_1.symlink_to(read_1)
                output_filename_2.symlink_to(read_2)
                logger.info(f"Original depth lower than {target_depth}. No downsampling needed.")

        return initial_estimated_vertical_depth

    @staticmethod
    def call_fastp(
        fastp_report: Path,
        input_read_1: Path,
        input_read_2: Path,
        output_read1: Path,
        output_read2: Path,
        logger: Logger,
    ):
        command_line = [
            f"{CONDA_BIN_DIR}/mamba", #Container
            "run", #Container
            "fastp",
            "--detect_adapter_for_pe",
        ]
        command_line += ["--in1", str(input_read_1)]
        command_line += ["--in2", str(input_read_2)]
        command_line += ["--out1", str(output_read1)]
        command_line += ["--out2", str(output_read2)]
        command_line += ["-j", str(fastp_report)]
        run_external(command_line, logger)

    @staticmethod
    def get_fastp_observations(fastp_report: Path) -> dict[str, float]:
        def as_percent(value: float | int) -> float:
            return round(value * 100, 2)

        def checkCycleKit(numCycles: int):

            kitList = {
                "300 cycle kit":150, 
                "500 cycle kit":250, 
                "600 cycle kit":300
            }
            cycle = "Non Standard"
            for k,v in kitList.items():
                if v-5 <= numCycles <= v+5:
                    cycle = k
                    break
            return cycle
            
        with open(fastp_report, encoding="utf-8") as reader:
            fastp_report_json = load(reader)
        observations: dict[str, float] = {}

        maxCyc = max(fastp_report_json['read1_before_filtering']['total_cycles'], fastp_report_json['read2_before_filtering']['total_cycles'])
        observations["inferredKit"] = checkCycleKit(maxCyc)

        observations["rawQ20"] = as_percent(
            fastp_report_json["summary"]["before_filtering"]["q20_rate"]
        )
        observations["rawQ30"] = as_percent(
            fastp_report_json["summary"]["before_filtering"]["q30_rate"]
        )
        observations["trimmedQ20"] = as_percent(
            fastp_report_json["summary"]["after_filtering"]["q20_rate"]
        )
        observations["trimmedQ30"] = as_percent(
            fastp_report_json["summary"]["after_filtering"]["q30_rate"]
        )
        observations["trimmedQ20R1"] = as_percent(
            fastp_report_json["read1_after_filtering"]["q20_bases"]
            / fastp_report_json["read1_after_filtering"]["total_bases"]
        )
        observations["trimmedQ20R2"] = as_percent(
            fastp_report_json["read2_after_filtering"]["q20_bases"]
            / fastp_report_json["read2_after_filtering"]["total_bases"]
        )
        observations["trimmedQ30R1"] = as_percent(
            fastp_report_json["read1_after_filtering"]["q30_bases"]
            / fastp_report_json["read1_after_filtering"]["total_bases"]
        )
        observations["trimmedQ30R2"] = as_percent(
            fastp_report_json["read2_after_filtering"]["q30_bases"]
            / fastp_report_json["read2_after_filtering"]["total_bases"]
        )
        observations["rawRead1AvgQC"] = round(sum(list(fastp_report_json["read1_before_filtering"]["quality_curves"]["mean"]))/len(list(fastp_report_json["read1_before_filtering"]["quality_curves"]["mean"])), 2)
        observations["rawRead2AvgQC"] = round(sum(list(fastp_report_json["read2_before_filtering"]["quality_curves"]["mean"]))/len(list(fastp_report_json["read2_before_filtering"]["quality_curves"]["mean"])), 2)
        observations["rawCombinedAvgQC"] = round((observations["rawRead1AvgQC"] + observations["rawRead2AvgQC"])/2, 2)
        
        observations["trimmedRead1AvgQC"] = round(sum(list(fastp_report_json["read1_after_filtering"]["quality_curves"]["mean"]))/len(list(fastp_report_json["read1_after_filtering"]["quality_curves"]["mean"])), 2)
        observations["trimmedRead2AvgQC"] = round(sum(list(fastp_report_json["read2_after_filtering"]["quality_curves"]["mean"]))/len(list(fastp_report_json["read2_after_filtering"]["quality_curves"]["mean"])), 2)
        observations["trimmedCombinedAvgQC"] = round((observations["trimmedRead1AvgQC"] + observations["trimmedRead2AvgQC"])/2, 2)
        
        observations["rawRead1MeanLength"] = fastp_report_json["summary"]["before_filtering"]["read1_mean_length"]
        observations["rawRead2MeanLength"] = fastp_report_json["summary"]["before_filtering"]["read2_mean_length"]
        observations["rawCombinedReadMeanLength"] = round((observations["rawRead1MeanLength"]+observations["rawRead2MeanLength"])/2, 2)

        observations["trimmedRead1MeanLength"] = fastp_report_json["summary"]["after_filtering"]["read1_mean_length"]
        observations["trimmedRead2MeanLength"] = fastp_report_json["summary"]["after_filtering"]["read2_mean_length"]
        observations["trimmedCombinedReadMeanLength"] = round((observations["trimmedRead1MeanLength"]+observations["trimmedRead2MeanLength"])/2, 2)

        return observations
