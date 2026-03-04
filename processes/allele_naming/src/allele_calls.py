from collections import defaultdict
from gzip import open as gzip_open
from json import load
from pathlib import Path
from xml.etree.ElementTree import parse

from src.models import AlleleCall, LocusCalls, LocusCallsList


class AlleleCallError(Exception):
    ...


def get_allele_calls_profile(allele_calls_profile_path: Path) -> list[str]:
    with gzip_open(allele_calls_profile_path, "rt", encoding="utf-8") as reader:
        allele_calls_profile = load(reader)
    # convert to list (dict.keys() returns an Iterable, which does not support all the features of a list)
    return list(allele_calls_profile["values"].keys())


def parse_allele_calls(allele_calls_file: Path) -> LocusCallsList:
    """
    Parse the allele calls XML file.
    Return a dictionary  with Locus ids and AlleleCalls

    Structure of the allele calls file:
    <loci>
        <l>
            <id>SALM_25362</id>                                         Locus identifier
            <as>                                                        A tag with one or more <a>
            <a>                                                         Contains all the info for a single allele call
                    <id val="7" />                                      Allele identifier
                    <si type="valid" val="100" />                       Sequence identity, equivalent to similarity
                    <ev type="valid" val="0" />                         E-value (from blast)
                    <bs type="valid" val="904" />                       Bitscore (from blast)
                    <no type="valid" val="0" />                         Number of “other” non-ACGT bases
                    <start type="info" val="99556" />                   Start location in the contig
                    <stop type="info" val="100057" />                   Stop location in the contig
                    <rs type="valid" val="0" />                         Repeat score: sequence identity of secondary hit
                    <al type="valid" val="501" />                       Alignment lengh (from blast)
                    <fwd type="info" val="1" />                         Fwd or rev
                    <cid type="info" val="Contig_34_43.018_pilon" />    Contig identifier
                    <nm type="valid" val="0" />                         Number of mismatches
                    <ngo type="valid" val="0" />                        Number of gaps open
                    <rss type="valid" val="N" />                        Requires start/stop
                    <bc type="valid" val="N" />                         Begin codon (=start codon)
                    <ec type="valid" val="N" />                         End codon
                    <fl type="valid" val="Y" />                         Full-length alignment
                    <is type="valid" val="N" />                         Internal stop
                    <s>ATTGC</s>                                        Sequence
                    <sb>CTGCGTGA</sb>                                   Before the sequence
                    <sa>GCCATCTC</sa>                                   After the sequence
                </a>
            </as>
        </l>
        ...
    </loci>

    The parser iterate over the different tags, <l>, then <as>, then <a> and finally all the tags in <a>.
    """
    with gzip_open(allele_calls_file, "rt", encoding="utf-8") as reader:
        root = parse(reader)

    locus_calls: LocusCalls = defaultdict(list)

    for locus in root.iter("l"):
        locus_id = locus.findtext("id")

        for allele in locus.iter("a"):
            allele_info = {
                node.tag: {
                    "val": node.attrib.get("val"),
                    "type": node.attrib.get("type"),
                    "text": node.text,
                }
                for node in allele
            }

            # an allele needs to have an id, a similarity and a sequence
            # probably it is legacy to check this, and they will always be present
            # for the BN CE assembly-free algorithm, the id can be "-2" to signal unknown
            # # but that is obsolete now
            # in fact, we will never use the id returned by the allele calling algorithm
            closest_allele_id = allele_info.get("id", {}).get("val")
            similarity = allele_info.get("si", {}).get("val")
            sequence = allele_info.get("s", {}).get("text")
            if not (closest_allele_id and similarity and sequence):
                continue

            sequences = {
                key: value["text"]
                for key, value in allele_info.items()
                if value["text"]
            }
            qualities = {
                key: value["val"]
                for key, value in allele_info.items()
                if value["type"] == "valid"
            }
            info = {
                key: value["val"]
                for key, value in allele_info.items()
                if value["type"] == "info"
            }

            locus_calls[locus_id].append(
                AlleleCall(
                    closest_allele_id=closest_allele_id,
                    sequences=sequences,
                    qualities=qualities,
                    info=info,
                )
            )

    return locus_calls


def filter_locus_calls(locus_calls: LocusCallsList, profile: list[str]) -> LocusCalls:
    """
    Filter locus calls to keep only allele calls that have been called by allele-calling
    """
    filtered_locus_calls: LocusCalls = {}

    for locus_id in profile:
        allele_calls = locus_calls.get(locus_id, [])

        # sanity check
        if len(allele_calls) != 1:
            raise AlleleCallError(
                f"Locus {locus_id} in profile has {len(allele_calls)} allele calls."
            )

        filtered_locus_calls[locus_id] = allele_calls[0]

    return filtered_locus_calls
