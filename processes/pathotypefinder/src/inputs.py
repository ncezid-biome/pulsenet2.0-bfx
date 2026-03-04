from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCInput,
    FastaInput,
    FastqInput,
    OrganismInput,
    Path    
)
from pydantic import Field


class PathotypeFinderInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    organism: OrganismInput = Field(
        description="Name of the organism genus."
    )
    reference: Path = Field(
        description="Path to reference fastas directory"
    )

