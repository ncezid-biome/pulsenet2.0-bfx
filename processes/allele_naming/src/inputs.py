from ngs_pipeline_lib.base.inputs import (
    BaseInputs,
    DirectoryPath,
    FilePath,
    JsonInput,
    KnowledgeBase,
)
from pydantic import Field


class AlleleCacheKB(KnowledgeBase):
    cache: DirectoryPath | None = Field(
        description="Path to the folder with the allele cache (accepted alleles fasta.gz)."
    )


class AlleleNamingInputs(BaseInputs):
    allele_calls_xml: FilePath = Field(
        description="Path to the BLASTAlleleFinder allele calls XML output file."
    )
    allele_calls_profile: FilePath = Field(
        description="Path to the allele calling profile standard_calls.json.gz file."
    )
    allele_cache_kb: AlleleCacheKB = Field(
        description="Path to the allele cache knowledge base folder."
    )
    nomenclature_settings: JsonInput = Field(
        description="Path to the nomenclature service settings json file."
    )
