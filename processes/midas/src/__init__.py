__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.midas import Midas
from src.inputs import MidasInputs


@cli.command(name="Midas")
def midas(inputs: MidasInputs):
    setup_logging(inputs.logging_dir)
    algorithm = Midas(inputs)
    algorithm.execute()
