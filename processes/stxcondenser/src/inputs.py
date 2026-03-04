from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCKnowledgeBase,
    FastqInput,
    DirectoryPath,
    FilePath
)
from pydantic import Field


class StxCondenserInputs(BaseInputs):
    stx_results: FilePath = Field(description="Path to stxfinder results file.")
    ispcr_results: FilePath = Field(description="Path to insilico PCR results file.")