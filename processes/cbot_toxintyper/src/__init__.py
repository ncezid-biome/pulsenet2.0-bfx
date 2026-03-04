__version__ = "1.0.0"

from ngs_pipeline_lib.cli import cli
from ngs_pipeline_lib.tools.logging import setup_logging

from src.cbot_toxintyping import CbotToxinTyping
from src.inputs import CbotToxinTypingInputs, ISPCRInputs
from src.ispcr import ISPCR

@cli.command(name="CbotToxinTyping")
def cbottoxintyping(inputs: CbotToxinTypingInputs):
    setup_logging(inputs.logging_dir)
    algorithm = CbotToxinTyping(inputs)
    algorithm.execute()

## Adding isPCR
@cli.command(name="isPCR")
def ispcr(inputs: ISPCRInputs):
    setup_logging(inputs.logging_dir)
    algorithm = ISPCR(inputs)
    algorithm.execute()