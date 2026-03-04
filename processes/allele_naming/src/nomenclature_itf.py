from dataclasses import dataclass
from xml.etree.ElementTree import Element, SubElement

from ngs_pipeline_lib.tools.tools import hash_sequence

from src.models import AlleleCall, LocusCalls, Profile
from src.nomenclature_service import NomenclatureService

YES_NO_MAP = {"N": "0", "Y": "1"}


def get_submit_allele_xml(locus_id: str, allele_call: AlleleCall) -> Element:
    """
    Get the XML format for allele submission
    """
    allele_node = Element("SubmitNewAllele")

    # main information
    SubElement(allele_node, "LocusId").text = locus_id
    SubElement(allele_node, "Id").text = allele_call.closest_allele_id
    SubElement(allele_node, "Hash").text = hash_sequence(allele_call.sequence)
    SubElement(allele_node, "Sequence").text = allele_call.sequence

    # additional information
    info_node = SubElement(allele_node, "Info")
    for key, value in allele_call.info.items():
        SubElement(info_node, key, type="i").text = value

    for key, value in allele_call.qualities.items():
        # quality (or "validation") values can be either Y|N or a float value
        # the Y|N needs to be mapped to 1|0 for allele submission
        # if the value is not a float (or int), then it can be skipped
        value = YES_NO_MAP.get(value, value)
        try:
            float(value)
        except ValueError:
            continue
        SubElement(info_node, key, type="v").text = value

    for key, value in allele_call.sequences.items():
        # do not add the actual allele sequence again
        if key == "s":
            continue
        SubElement(info_node, key).text = value

    return allele_node


def get_submit_alleles_xml(
    organism_id: str, lab_id: str, sample_id: str, locus_calls: LocusCalls
) -> Element:
    """
    Get the XML format for alleles submission
    """
    submit_alleles = Element("SubmitNewAlleles")
    SubElement(submit_alleles, "Organism").text = organism_id
    SubElement(submit_alleles, "LabId").text = lab_id
    SubElement(submit_alleles, "Entry").text = sample_id
    SubElement(
        submit_alleles, "Comment"
    ).text = f"New alleles for {sample_id} submitted by {lab_id}"
    submit_alleles.extend(
        get_submit_allele_xml(locus_id, allele_call)
        for locus_id, allele_call in locus_calls.items()
    )
    return submit_alleles


def iter_alleles(locus_calls: LocusCalls) -> LocusCalls:
    """
    Iterate over alleles to submit and return in batches of 50.
    Submitting alleles in batches of 50 gives the best performance;
    the best balance between making many calls (submitting one-by-one)
    and long processing on the server (submitting all at once).
    """
    BATCH_SIZE = 50
    batch = LocusCalls()
    for locus_id, allele_call in locus_calls.items():
        batch[locus_id] = allele_call
        if len(batch) == BATCH_SIZE:
            yield batch
            batch = LocusCalls()
    if batch:
        yield batch


@dataclass
class NomenclatureInterface:
    nomenclature_service: NomenclatureService
    organism_id: str
    lab_id: str
    sample_id: str

    def submit_alleles(self, alleles_to_submit: LocusCalls) -> Profile:
        """
        Submit alleles to nomenclature service to get allele ids
        Quality filtering has already been done before
        There is only one allele per locus
        """
        # sanity check
        organism_ids = self.nomenclature_service.get_organisms()
        if self.organism_id not in organism_ids:
            raise RuntimeError(
                f"Invalid nomenclature configuration. Organism {self.organism_id} not in available organism ids {organism_ids}."
            )

        profile = Profile()
        for allele_batch in iter_alleles(alleles_to_submit):
            named_alleles = self.nomenclature_service.submit_alleles(
                get_submit_alleles_xml(
                    self.organism_id, self.lab_id, self.sample_id, allele_batch
                )
            )
            profile.update(
                {
                    allele.findtext("LocusId"): allele.findtext("AlleleId")
                    for allele in named_alleles.iter("ResponseNewAllele")
                }
            )

        return profile
