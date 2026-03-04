__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.stx_condenser import StxCondenser
from src.inputs import StxCondenserInputs


@cli.command(name="stx_condenser")
def stxcondenser(inputs: StxCondenserInputs):
    setup_logging(inputs.logging_dir)
    algorithm = StxCondenser(inputs)
    algorithm.execute()
