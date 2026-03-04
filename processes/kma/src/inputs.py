from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCKnowledgeBase,
    FastqInput,
    FilePath
)
from pydantic import Field


class EcKMASerotyperInputs(BaseInputs):
    read1: FastqInput = Field(description="Filepath to the first read file.")
    read2: FastqInput = Field(description="Filepath to the second read file.")
    otypes: FilePath = Field(description="O type fsa reference file.")
    htypes: FilePath = Field(description="H type fsa reference file.")
    lookuptable: FilePath = Field(description="Filepath to custom O-group translation table.")
