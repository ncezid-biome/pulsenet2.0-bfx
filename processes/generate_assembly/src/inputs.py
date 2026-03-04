from ngs_pipeline_lib.base.inputs import BaseInputs, FastqInput, OrganismInput, QCKnowledgeBase
from pydantic import Field


class GenerateAssemblyInputs(BaseInputs):
    read1: FastqInput = Field(description="Path to the first read file.")
    read2: FastqInput = Field(description="Path to the second read file.")
    organism: OrganismInput = Field(
        description="Organism, either genus or genus and species."
    )
    qc_kb: QCKnowledgeBase = Field(description="Path to quality control knowledge base folder.")
    
