from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    FastaInput,
    JsonInput,
    KnowledgeBase,
    OrganismInput,
    QCKnowledgeBase,
)
from pydantic import Field


class SchemeMappingKB(KnowledgeBase):
    mapping: JsonInput | None = Field(
        description="Path to the file with organism -> scheme mapping."
    )


class MLSTInputs(BaseInputs):
    assembly: FastaInput = Field(description="Filepath to the assembly file.")
    qc_kb: QCKnowledgeBase = Field(
        description="Path to the quality control knowledge base folder."
    )
    scheme_mapping_kb: SchemeMappingKB = Field(
        description="Path to the Organism->Scheme mapping KB"
    )
    organism: OrganismInput = Field(
        description="Name of the expected organism, genus + (species).",
        default=OrganismInput(genus=None, species=None),
    )
