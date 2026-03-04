__version__ = "1.1.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.seqsero import SeqSero
from src.inputs import SeqSeroInputs


@cli.command(name="seqsero")
def seqsero(inputs: SeqSeroInputs):
    setup_logging(inputs.logging_dir)
    algorithm = SeqSero(inputs)
    algorithm.execute()

