from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCKnowledgeBase,
    FastqInput,
    DirectoryPath
)
from pydantic import Field


class EcStxTyperInputs(BaseInputs):
    read1: FastqInput = Field(description="Filepath to the first read file.")
    read2: FastqInput = Field(description="Filepath to the second read file.")
    holotoxins: DirectoryPath = Field(description="DirectoryPath to stx holotoxin fsa reference file(s).")