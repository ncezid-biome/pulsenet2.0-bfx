__version__ = "3.0.2"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.correct_assembly import CorrectAssembly
from src.inputs import CorrectAssemblyInputs


@cli.command(name="CorrectAssembly")
def correct_assembly(args: CorrectAssemblyInputs):
    setup_logging(args.logging_dir)
    algorithm = CorrectAssembly(args)
    algorithm.execute()
