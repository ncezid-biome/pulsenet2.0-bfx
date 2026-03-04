__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.plasmidfinder import PlasmidFinder
from src.inputs import PlasmidFinderInputs


@cli.command(name="PlasmidFinder")
def plasmidfinder(inputs: PlasmidFinderInputs):
    setup_logging(inputs.logging_dir)
    algorithm = PlasmidFinder(inputs)
    algorithm.execute()

