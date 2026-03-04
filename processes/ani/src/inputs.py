from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCInput,
    FastaInput,
    FilePath,
    
)
from pydantic import Field


class ANIInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    reference: FilePath = Field(
        description="Path to references.tsv"
    )
    qc: QCInput = Field(
        description="Path to quality control reference file."
    )
    qc_section: str = Field(
        description="QC rules section to use. 'ANI' or 'ORG_ANI'"
    )

