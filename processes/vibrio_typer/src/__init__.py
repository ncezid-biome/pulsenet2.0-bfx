__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.vibrio_virulence_finder import VibrioVirulenceFinder
from src.inputs import VibrioVirulenceInputs


@cli.command(name="vibrio_virulence")
def vibrio_virulence(inputs: VibrioVirulenceInputs):
    setup_logging(inputs.logging_dir)
    algorithm = VibrioVirulenceFinder(inputs)
    algorithm.execute()
