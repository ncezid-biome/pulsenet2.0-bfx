from dataclasses import dataclass


class Profile(dict[str, str]):
    ...


@dataclass
class AlleleCall:
    closest_allele_id: str
    qualities: dict[str, str]
    info: dict[str, str]
    sequences: dict[str, str]

    @property
    def allele_id(self) -> str | None:
        """
        Only return an allele id if the closest found allele is an exact match
        """
        return self.closest_allele_id if self.similarity == 100.0 else None

    @property
    def sequence(self) -> str:
        return self.sequences["s"]

    @property
    def similarity(self) -> float:
        return float(self.qualities.get("si", 0.0))


class LocusCallsList(dict[str, list[AlleleCall]]):
    ...


class LocusCalls(dict[str, AlleleCall]):
    ...
