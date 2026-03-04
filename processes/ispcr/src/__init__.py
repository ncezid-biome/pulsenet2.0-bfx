__version__ = "1.0.2"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.ispcr import ISPCR
from src.inputs import ISPCRInputs


@cli.command(name="isPCR")
def ispcr(inputs: ISPCRInputs):
    setup_logging(inputs.logging_dir)
    algorithm = ISPCR(inputs)
    algorithm.execute()

