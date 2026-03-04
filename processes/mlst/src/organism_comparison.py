from enum import IntEnum
from os import getenv

from ngs_pipeline_lib.base.inputs import OrganismInput

MLST_SCHEME_SPECIES_MAP_PATH = getenv(
    "MLST_SCHEME_SPECIES_MAP_PATH", "/opt/conda/db/scheme_species_map.tab"
)
#MLST_SCHEME_SPECIES_MAP_PATH = "/scicomp/home-pure/lsz0/PN2.0/NewStructure/MLST/pn2.0_test/scheme_species_map.tab" ## no container


class ComparisonResult(IntEnum):
    ok = 0
    different_species = 1
    different_genus = 2
    observed_organism_unknown = 3

    def __str__(self):
        return str(self.name)


def compare_organism(
    observed_organism: OrganismInput | None,
    process_organism: OrganismInput,
) -> ComparisonResult:
    """
    process_organism and observed_organism always have at least a genus
    """
    if observed_organism is None:
        return ComparisonResult.observed_organism_unknown
    if observed_organism.genus != process_organism.genus:
        return ComparisonResult.different_genus
    if (
        observed_organism.species is not None
        and process_organism.species is not None
        and observed_organism.species != process_organism.species
    ):
        return ComparisonResult.different_species
    return ComparisonResult.ok
