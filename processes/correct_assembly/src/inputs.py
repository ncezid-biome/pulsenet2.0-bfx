from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    FastaInput,
    FastqInput,
    JsonInput,
    KnowledgeBase,
    OrganismInput,
    QCKnowledgeBase,
)
from pydantic import Field


class CleaningKnowledgeBase(KnowledgeBase):
    settings: JsonInput | None = Field(
        description="Path to the organism-to-cleaning-settings mapping file."
    )


class CorrectAssemblyInputs(BaseInputs):
    organism: OrganismInput = Field(
        description="Organism, either genus or genus and species."
    )
    assembly: FastaInput = Field(description="Path to the assembly file.")
    read1: FastqInput = Field(description="Path to the first read file.")
    read2: FastqInput = Field(description="Path to the second read file.")
    cleaning_kb: CleaningKnowledgeBase = Field(
        description="Path to cleaning knowledge base folder.",
    )
    qc_kb: QCKnowledgeBase = Field(
        description="Path to quality control knowledge base folder."
    )
    bowtie: bool = Field(description="Whether to use bowtie (or bwa).", default=False)
    bowtie_n_hits: int = Field(description="Number of hits for bowtie.", default=1)
