__version__ = "3.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.generate_assembly import GenerateAssembly
from src.inputs import GenerateAssemblyInputs


@cli.command(name="GenerateAssembly")
def generate_assembly(args: GenerateAssemblyInputs):
    setup_logging(args.logging_dir)
    algorithm = GenerateAssembly(args)
    algorithm.execute()
