__version__ = "3.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.cleanup_reads import CleanupReads
from src.download_reads import DownloadReads
from src.inputs import CleanupReadsInputs, DownloadReadsInputs


@cli.command(name="DownloadReads")
def download_reads(args: DownloadReadsInputs):
    setup_logging(args.logging_dir)
    algorithm = DownloadReads(args)
    algorithm.execute()


@cli.command(name="CleanupReads")
def cleanup_reads(args: CleanupReadsInputs):
    setup_logging(args.logging_dir)
    algorithm = CleanupReads(args)
    algorithm.execute()
