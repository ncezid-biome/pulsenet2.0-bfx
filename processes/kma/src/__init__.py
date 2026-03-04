__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.ec_kma_serotyper import EcKMASerotyper
from src.inputs import EcKMASerotyperInputs


@cli.command(name="ec_kma_serotyper")
def ec_kma_serotyper(inputs: EcKMASerotyperInputs):
    setup_logging(inputs.logging_dir)
    algorithm = EcKMASerotyper(inputs)
    algorithm.execute()
