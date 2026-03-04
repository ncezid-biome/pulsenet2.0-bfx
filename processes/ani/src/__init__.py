__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.ani import ANI
from src.inputs import ANIInputs


@cli.command(name="ANI")
def ani(inputs: ANIInputs):
    setup_logging(inputs.logging_dir)
    algorithm = ANI(inputs)
    algorithm.execute()

