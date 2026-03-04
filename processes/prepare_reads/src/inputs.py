from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    FastqInput,
    JsonInput,
    KnowledgeBase,
    OrganismInput,
    QCKnowledgeBase,
)
from pydantic import Field


class GenomeKnowledgeBase(KnowledgeBase):
    sizes: JsonInput | None = Field(
        description="Path to the organism-to-genome-size mapping file."
    )
    target_depths: JsonInput | None = Field(
        description="Path to the organism-to-target-depth mapping file."
    )


class DownloadReadsInputs(BaseInputs):
    aws: bool = Field(description="Whether or not to first try to download from AWS.")
    accession_id: str = Field(description="NCBI accession id.")


class CleanupReadsInputs(BaseInputs):
    read1: FastqInput = Field(description="Path to the first read file.")
    read2: FastqInput = Field(description="Path to the second read file.")
    organism: OrganismInput = Field(
        description="Name of the organism, either genus or genus + species."
    )
    genome_kb: GenomeKnowledgeBase = Field(
        description="Path to the genome knowledge base folder."
    )
    qc_kb: QCKnowledgeBase = Field(
        description="Path to the quality control knowledge base folder."
    )
