__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.ec_stxtyper import EcStxTyper
from src.inputs import EcStxTyperInputs


@cli.command(name="ec_stxtyper")
def ec_stxtyper(inputs: EcStxTyperInputs):
    setup_logging(inputs.logging_dir)
    algorithm = EcStxTyper(inputs)
    algorithm.execute()
