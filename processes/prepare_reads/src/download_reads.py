import os
from logging import Logger
from os import environ
from pathlib import Path
from shutil import rmtree

from boto3 import resource
from boto3.exceptions import RetriesExceededError
from botocore.exceptions import ClientError
from botocore.handlers import disable_signing
from mypy_boto3_s3 import S3ServiceResource
from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.tools.stub import generate_fake_fastq_file

from src.inputs import DownloadReadsInputs
from src.outputs import DownloadReadsOutputs

_SRA_BUCKET = "sra-pub-run-odp"
CONDA_BIN_DIR = os.getenv("CONDA_BIN_DIR", "/opt/conda/condabin")


class DownloadReads(Algorithm[DownloadReadsInputs, DownloadReadsOutputs]):
    """
    See downloadFastq from reads.nf
    """

    outputs_class = DownloadReadsOutputs

    def execute_stub(self):
        self.outputs.read_1.content = generate_fake_fastq_file()
        self.outputs.read_2.content = generate_fake_fastq_file()

        self.set_results()

    def execute_implementation(self) -> None:
        self.logger.info(environ["HOME"])
        # uses prefetch if downloading from AWS is disabled or failed
        if not (
            self.inputs.aws
            and self.download_from_aws(self.inputs.accession_id, self.logger)
        ):
            self.prefetch(self.inputs.accession_id, self.logger)

        self.split_sra_file(self.inputs.accession_id, self.logger)
        self.rename_fastq_files(self.inputs.accession_id)
        self.cleanup(self.inputs.accession_id)

        self.set_results()

    def set_results(self):
        self.result = {
            "read_1": self.inputs.publish_dir + str(self.outputs.read_1.path),
            "read_2": self.inputs.publish_dir + str(self.outputs.read_2.path),
        }

    @staticmethod
    def download_from_aws(accession_id: str, logger: Logger) -> bool:
        try:
            s3_resource: S3ServiceResource = resource("s3")
            s3_resource.meta.client.meta.events.register(
                "choose-signer.s3.*", disable_signing
            )
            key = f"sra/{accession_id}/{accession_id}"
            s3_resource.Bucket(_SRA_BUCKET).download_file(key, str(accession_id))
        except (RetriesExceededError, ClientError) as error:
            logger.error(error)
            return False
        return True

    @staticmethod
    def prefetch(accession_id: str, logger: Logger):
        run_external(
            [f"{CONDA_BIN_DIR}/mamba", "run", "prefetch", accession_id, "-O", "."],
            logger,
        )

    @staticmethod
    def split_sra_file(accession_id: str, logger: Logger):
        run_external(
            [
                f"{CONDA_BIN_DIR}/mamba",
                "run",
                "fasterq-dump",
                "--split-3",
                accession_id,
            ],
            logger,
        )

    @staticmethod
    def rename_fastq_files(accession_id: str):
        suffixes = ["_1", "_2"]
        for suffix in suffixes:
            Path(f"{accession_id}{suffix}.fastq").rename(f"read{suffix}.fastq")

    @staticmethod
    def cleanup(accession_id: str):
        rmtree(Path() / ".ncbi", ignore_errors=True)
        rmtree(accession_id, ignore_errors=True)
