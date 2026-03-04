__version__ = "1.1.2"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.inputs import MLSTInputs
from src.mlst import MLST


@cli.command(name="MLST")
def mlst_command(args: MLSTInputs):
    setup_logging(args.logging_dir)
    algorithm = MLST(args)
    algorithm.execute()
