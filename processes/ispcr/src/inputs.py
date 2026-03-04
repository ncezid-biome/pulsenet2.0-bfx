from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCInput,
    FastaInput,
    FastqInput,
    OrganismInput,
    FilePath,
    
)
from pydantic import Field


class ISPCRInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    organism: OrganismInput = Field(
        description="Name of the organism, either genus or genus + species."
    )

