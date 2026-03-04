__version__ = "3.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.allele_naming import AlleleNaming
from src.inputs import AlleleNamingInputs


@cli.command(name="AlleleNaming")
def allele_naming(inputs: AlleleNamingInputs):
    setup_logging(inputs.logging_dir)
    algorithm = AlleleNaming(inputs)
    algorithm.execute()
