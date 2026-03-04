from enum import IntEnum
from typing import Any, Self

from src.organism_comparison import ComparisonResult

MISSING_ALLELES_FOR_ERROR = 2
TOTAL_ALLELE_ISSUES_FOR_ERROR = 5
MISSING_ALLELES_FOR_WARNING = 1
MULTIPLE_ALLELES_FOR_WARNING = 2
TOTAL_ALLELE_ISSUES_FOR_WARNING = 2


class AlleleMatch(IntEnum):
    exact = 0
    novel = 1
    partial = 2
    missing = 3
    unknown = 4

    def __str__(self):
        return str(self.name)


def check_allele(allele: str) -> list[AlleleMatch]:
    def _check_allele(allele: str) -> AlleleMatch:
        if allele.isdigit():
            return AlleleMatch.exact
        if allele.startswith("~") and allele[1:].isdigit():
            return AlleleMatch.novel
        if allele.endswith("?") and allele[:1].isdigit():
            return AlleleMatch.partial
        if allele == "-":
            return AlleleMatch.missing
        return AlleleMatch.unknown

    return [
        _check_allele(each_allele.strip()) for each_allele in allele.strip().split(",")
    ]


# missing: 1 -> warning *
# missing 2 or more -> fail
# multiple: 2 or more -> warning (possible contamination)
# 1 problem -> ok
# 2 to 4 problems -> warning *
# 5 problems -> fail


class AlleleIssues:
    def __init__(
        self,
        scheme: str | None = None,
        sequence_type: int | None = None,
        genes: list[str] | None = None,
        exact_alleles: list[str] | None = None,
        unknown_alleles: list[str] | None = None,
        missing_alleles: list[str] | None = None,
        multiple_alleles: list[str] | None = None,
        novel_alleles: list[str] | None = None,
        partial_alleles: list[str] | None = None,
        organism_comparison: ComparisonResult | None = None,
    ):
        self.scheme = scheme if scheme != "-" else None
        self.sequence_type: int | None = sequence_type
        self.genes: list[str] = genes or []

        self.organism_comparison: ComparisonResult = (
            organism_comparison
            if organism_comparison is not None
            else ComparisonResult.observed_organism_unknown
        )

        self.exact_alleles: list[str] = exact_alleles or []
        self.unknown_alleles: list[str] = unknown_alleles or []
        self.missing_alleles: list[str] = missing_alleles or []
        self.multiple_alleles: list[str] = multiple_alleles or []
        self.novel_alleles: list[str] = novel_alleles or []
        self.partial_alleles: list[str] = partial_alleles or []

    @classmethod
    def from_scheme_result(
        cls, scheme_result: dict[str, Any], organism_comparison: ComparisonResult | None
    ) -> Self:
        scheme = (
            scheme
            if ((scheme := scheme_result.get("scheme")) and scheme != "-")
            else None
        )
        sequence_type = (
            int(sequence_type)
            if (
                (sequence_type := scheme_result.get("sequence_type"))
                and sequence_type != "-"
            )
            else None
        )

        gene_to_alleles = {
            gene: check_allele(allele)
            for (gene, allele) in scheme_result.get("alleles", {}).items()
        }

        exact_alleles = [
            gene
            for (gene, alleles) in gene_to_alleles.items()
            if max(alleles) == AlleleMatch.exact
        ]
        unknown_alleles = [
            gene
            for (gene, alleles) in gene_to_alleles.items()
            if max(alleles) == AlleleMatch.unknown
        ]
        missing_alleles = [
            gene
            for (gene, alleles) in gene_to_alleles.items()
            if max(alleles) == AlleleMatch.missing
        ]
        multiple_alleles = [
            gene for (gene, alleles) in gene_to_alleles.items() if len(alleles) > 1
        ]
        novel_alleles = [
            gene
            for (gene, alleles) in gene_to_alleles.items()
            if max(alleles) == AlleleMatch.novel
        ]
        partial_alleles = [
            gene
            for (gene, alleles) in gene_to_alleles.items()
            if max(alleles) == AlleleMatch.partial
        ]

        return cls(
            scheme=scheme,
            sequence_type=sequence_type,
            genes=list(gene_to_alleles.keys()),
            exact_alleles=exact_alleles,
            unknown_alleles=unknown_alleles,
            missing_alleles=missing_alleles,
            multiple_alleles=multiple_alleles,
            novel_alleles=novel_alleles,
            partial_alleles=partial_alleles,
            organism_comparison=organism_comparison,
        )

    @property
    def total_issues(self) -> int:
        return len(
            self.unknown_alleles
            + self.missing_alleles
            + self.multiple_alleles
            + self.novel_alleles
            + self.partial_alleles
        )

    @property
    def all_issues(self) -> int:
        return list(
            self.unknown_alleles
            + self.missing_alleles
            + self.multiple_alleles
            + self.novel_alleles
            + self.partial_alleles
        )

    @property
    def has_error(self) -> bool:
        if not self.sequence_type and (
            not self.genes
            or len(self.unknown_alleles) > 0
            or len(self.missing_alleles) >= MISSING_ALLELES_FOR_ERROR
            or self.total_issues >= TOTAL_ALLELE_ISSUES_FOR_ERROR
        ):
            return True
        return False

    @property
    def has_warning(self) -> bool:
        if self.has_error or (
            not self.sequence_type
            and (
                not self.genes
                or len(self.missing_alleles) >= MISSING_ALLELES_FOR_WARNING
                or len(self.multiple_alleles) >= MULTIPLE_ALLELES_FOR_WARNING
                or self.total_issues >= TOTAL_ALLELE_ISSUES_FOR_WARNING
            )
        ):
            return True
        return False

    def __repr__(self) -> str:
        return (
            "<AlleleIssues "
            f"scheme:{self.scheme}, "
            f"ST:{self.sequence_type}, "
            f"Warn:{self.has_warning}, "
            f"Error:{self.has_error}, "
            f"Issues:{self.total_issues}, "
            ">"
        )

    def __lt__(self, other: Self) -> bool:
        """
        Considers the least errors/problems as the 'lower' one, so that when doing:
            min([AlleleIssues1, AlleleIssues2, AlleleIssues3])
        this will return the one with the least amount of issues
        """
        # if one has ST and the other does not, the ST is taken as 'lower'
        if self.sequence_type is not None and other.sequence_type is None:
            return True
        if other.sequence_type is not None and self.sequence_type is None:
            return False

        # if both have STs, the closest (organism) ComparisonResult is 'lower'
        if self.sequence_type is not None and other.sequence_type is not None:
            if (
                self.organism_comparison != other.organism_comparison
                and ComparisonResult.observed_organism_unknown
                not in [self.organism_comparison, other.organism_comparison]
            ):
                return self.organism_comparison < other.organism_comparison

        # If both don't have STs or both have and have the same organism in comparison,
        # checks for errors, and will consider the one without errors as 'lower'
        # than one with errors
        if self.has_error and not other.has_error:
            return False
        if other.has_error and not self.has_error:
            return True

        # Checks for warnings, and will consider the one without warnings as
        # 'lower' than one with warnings
        if self.has_warning and not other.has_warning:
            return False
        if other.has_warning and not self.has_warning:
            return True

        # Checks for unknown alleles, and will consider the one with the least
        # unknown alleles as 'lower'
        if len(self.unknown_alleles) > len(other.unknown_alleles):
            return False
        if len(other.unknown_alleles) > len(self.unknown_alleles):
            return True

        # Checks for missing alleles, and will consider the one with the least
        # missing alleles as 'lower'
        if len(self.missing_alleles) > len(other.missing_alleles):
            return False
        if len(other.missing_alleles) > len(self.missing_alleles):
            return True

        # Checks for the organism comparison, and will consider the one with the
        # closest comparison as 'lower'
        if self.organism_comparison != other.organism_comparison:
            return self.organism_comparison < other.organism_comparison

        # Checks for total issues, and will consider the one with the least
        # total issues as 'lower'
        if self.total_issues > other.total_issues:
            return False
        if other.total_issues > self.total_issues:
            return True

        # Checks again for total issues, but this time ignoring novel alleles,
        # which are a minor issue, and will consider the one with the least
        # issues as 'lower'
        if (self.total_issues - len(self.novel_alleles)) > (
            other.total_issues - len(other.novel_alleles)
        ):
            return False
        if (other.total_issues - len(other.novel_alleles)) > (
            self.total_issues - len(self.novel_alleles)
        ):
            return True

        # If everything else fails, tries to return the one with the most
        # exact alleles as 'lower'. Notice that if both have the same number,
        # it will return False
        return len(self.exact_alleles) > len(other.exact_alleles)

    # def __eq__(self, other) -> bool:
    #     return (
    #         str(self) == str(other)
    #         and self.sequence_type == other.sequence_type
    #         and self.genes == other.genes
    #         and self.has_error == other.has_error
    #         and self.has_warning == other.has_warning
    #         and self.total_issues == other.total_issues
    #         and self.exact_alleles == other.exact_alleles
    #         and self.missing_alleles == other.missing_alleles
    #         and self.multiple_alleles == other.multiple_alleles
    #         and self.novel_alleles == other.novel_alleles
    #         and self.partial_alleles == other.partial_alleles
    #         and self.organism_comparison == other.organism_comparison
    #     )

    def __str__(self) -> str:
        return (
            "<Allele issues ("
            f"ST: {self.sequence_type}, "
            f"organism_comparison: {self.organism_comparison}, "
            f"missing: {len(self.missing_alleles)}, "
            f"multiple: {len(self.missing_alleles)}, "
            f"total_issues: {self.total_issues}"
            ")>"
        )
