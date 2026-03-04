from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.base.outputs import BaseOutputs

from src.allele_cache import AlleleCache
from src.allele_calls import (
    filter_locus_calls,
    get_allele_calls_profile,
    parse_allele_calls,
)
from src.inputs import AlleleNamingInputs
from src.models import LocusCalls, Profile
from src.nomenclature_itf import NomenclatureInterface
from src.nomenclature_service import NomenclatureService


class AlleleNaming(Algorithm[AlleleNamingInputs, BaseOutputs]):
    outputs_class = BaseOutputs

    def execute_stub(self):
        """
        Execute in stub mode.
        """
        allele_calling_profile = get_allele_calls_profile(
            self.inputs.allele_calls_profile
        )
        allele_naming_profile = {locus_id: "1" for locus_id in allele_calling_profile}
        self.result = {"values": allele_naming_profile}
        self.qc_report.add_metric("submitted_alleles_count", 0)

    def execute_implementation(self) -> None:
        """
        Parses the allele_calls_output and looks up the called allele sequences in allele_calls_xml.
        The allele sequences are searched in the known alleles,
        and unknown alleles are submitted to the nomenclature service.
        """

        allele_calling_profile = get_allele_calls_profile(
            self.inputs.allele_calls_profile
        )
        allele_calling_details = parse_allele_calls(self.inputs.allele_calls_xml)
        allele_calls = filter_locus_calls(
            allele_calling_details, allele_calling_profile
        )

        wgmlst_profile, alleles_to_submit = self.get_allele_names(allele_calls)
        if alleles_to_submit:
            submitted_alleles = self.submit_alleles(alleles_to_submit)
            wgmlst_profile.update(submitted_alleles)
        self.qc_report.add_metric("submitted_alleles_count", len(alleles_to_submit))

        self.result = {"values": wgmlst_profile}

    def get_allele_names(self, locus_calls: LocusCalls) -> tuple[Profile, LocusCalls]:
        """
        Lookup allele names in the cache
        """
        allele_cache = AlleleCache(self.inputs.allele_cache_kb.cache)
        known_alleles = Profile()
        alleles_to_submit = LocusCalls()
        for locus_id, allele_call in locus_calls.items():
            if allele_call.allele_id:
                allele_id = allele_call.allele_id
            else:
                allele_id = allele_cache.get_allele_id(locus_id, allele_call.sequence)
            if allele_id:
                known_alleles[locus_id] = allele_id
            else:
                alleles_to_submit[locus_id] = allele_call
        return known_alleles, alleles_to_submit

    def submit_alleles(self, alleles_to_submit: LocusCalls) -> Profile:
        """
        Submit new alleles to nomenclature service
        """
        settings = self.inputs.nomenclature_settings.get_dict()
        nomenclature_service = NomenclatureService(
            settings["url"],
            settings["project"],
            settings["password"],
            settings["serial"],
        )
        nomenclature_itf = NomenclatureInterface(
            nomenclature_service,
            settings["organism_id"],
            settings.get("lab_id", "lab_id"),
            self.inputs.sample_id,
        )
        return nomenclature_itf.submit_alleles(alleles_to_submit)
