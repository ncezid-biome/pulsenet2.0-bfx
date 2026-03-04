from ngs_pipeline_lib.base.inputs import BaseInputs, FastaInput, FilePath, KnowledgeBase, OrganismInput
from pydantic import Field


class AMRFinderKB(KnowledgeBase):
    db: FilePath | None = Field(
        description="Path to the AMRFinder database tar.gz file."
    )


class FindGenesInputs(BaseInputs):
    assembly: FastaInput = Field(description="Path to the assembly file.")
    organism: OrganismInput = Field(
        description="Name of the organism, either genus or genus + species."
    )

