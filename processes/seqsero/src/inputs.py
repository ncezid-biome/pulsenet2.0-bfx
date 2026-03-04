from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    QCInput,
    FastaInput,
    FastqInput,
    OrganismInput,
    FilePath,
    
)
from pydantic import Field

class SeqSeroInputs(BaseInputs):
    assembly: FastaInput = Field(
        description="Path to assembly file."
    )
    organism: OrganismInput = Field(
        description="Name of the organism, either genus or genus + species."
    )
    read1: FastqInput = Field(
        description="Filepath to the first read file."
    )
    read2: FastqInput = Field(
        description="Filepath to the second read file."
    )
    ispcr_out: FilePath = Field(
        description="Path to isPCR output file"
    )
    species_ani_out: FilePath = Field(
        description="Path to species ANI best_hit.json output file"
    )
    mlst_out: FilePath = Field(
        description="Path to MLST outputs.json."
    )

