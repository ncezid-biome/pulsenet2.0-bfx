__version__ = "3.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.find_genes import FindGenes
from src.inputs import FindGenesInputs


@cli.command(name="FindGenes")
def find_genes(args: FindGenesInputs):
    setup_logging(args.logging_dir)
    algorithm = FindGenes(args)
    algorithm.execute()
