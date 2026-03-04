__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.pathotypefinder import PathotypeFinder
from src.inputs import PathotypeFinderInputs


@cli.command(name="PathotypeFinder")
def pathotypefinder(inputs: PathotypeFinderInputs):
    setup_logging(inputs.logging_dir)
    algorithm = PathotypeFinder(inputs)
    algorithm.execute()
