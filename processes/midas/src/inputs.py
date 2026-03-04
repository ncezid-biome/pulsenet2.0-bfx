from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCKnowledgeBase,
    FastqInput,
    OrganismInput,
)
from pydantic import Field


class MidasInputs(BaseInputs):
    read1: FastqInput = Field(description="Filepath to the first read file.")
    read2: FastqInput = Field(description="Filepath to the second read file.")
    qc_kb: QCKnowledgeBase = Field(
        description="Path to quality control knowledge base folder."
    )
    organism: OrganismInput = Field(
        description="Organism, either genus or genus and species."
    )