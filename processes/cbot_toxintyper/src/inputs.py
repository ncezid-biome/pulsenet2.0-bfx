from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCInput,
    FastaInput,
    FastqInput,
    OrganismInput,
    Path    
)
from pydantic import Field


class CbotToxinTypingInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    organism: OrganismInput = Field(
        description="Name of the organism genus."
    )
    reference: Path = Field(
        description="Path to reference fastas directory"
    )

class ISPCRInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    organism: OrganismInput = Field(
        description="Name of the organism, either genus or genus + species."
    )