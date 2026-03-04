from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCKnowledgeBase,
    FastaInput,
    FastqInput,
    OrganismInput,
    DirectoryPath,
    FilePath
)
from pydantic import Field


class VibrioVirulenceInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    organism: OrganismInput = Field(
        description="Name of the organism, either genus or genus + species"
    )
