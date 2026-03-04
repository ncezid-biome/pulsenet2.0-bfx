from collections import defaultdict
from gzip import open as gzip_open
from pathlib import Path

from Bio.SeqIO import parse


class AlleleCacheError(Exception):
    ...


class AlleleCache:
    def __init__(self, allele_cache_dir: Path):
        """
        Load the accepted alleles cache
        Beware that an allele sequence might appear in more than one locus (should not, but you never know...)
        """
        self._cache: dict[str, dict[str, str]] = defaultdict(dict)

        accepted_alleles_link_file = allele_cache_dir / "acceptedalleles_link"
        if not accepted_alleles_link_file.is_file():
            raise AlleleCacheError(
                "Invalid allele cache: no accepted alleles link file found."
            )

        with open(
            accepted_alleles_link_file, encoding="utf-8"
        ) as accepted_alleles_link_reader:
            accepted_alleles_link = accepted_alleles_link_reader.read(1)

        allele_cache_file = (
            allele_cache_dir / f"acceptedalleles_{accepted_alleles_link}.fasta.gz"
        )
        if not allele_cache_file.is_file():
            raise AlleleCacheError(
                "Invalid allele cache: the accepted alleles fasta file is not found."
            )

        with gzip_open(allele_cache_file, "rt", encoding="utf-8") as allele_reader:
            for allele in parse(allele_reader, "fasta"):
                # the allele.name is constructed like "<locus_id>_<allele_id>"
                # the allele_id is an integer number, and the locus_id is constructed like "<organism_abbr>_<locus_nr>"
                # E.g. "LMO_1_1" -> locus_id = "LMO_1" and allele_id = "1"
                locus_id, allele_id = allele.name.rsplit("_", 1)
                # make the sequence upper case for case-insensitive comparison
                allele_seq = str(allele.seq).upper()
                self._cache[locus_id][allele_seq] = allele_id

    def get_allele_id(self, locus_id: str, allele_sequence: str) -> str:
        """
        Get the allele id for an allele sequence of a locus
        Return None if the locus or allele is not found
        """
        return self._cache.get(locus_id, {}).get(allele_sequence.upper())
