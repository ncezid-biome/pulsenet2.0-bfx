from csv import reader as csv_reader
from json import load as json_load
from logging import Logger
from os import getenv, makedirs
from pathlib import Path
from typing import Iterable

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.base.inputs import OrganismInput
from ngs_pipeline_lib.base.report import QCReport
from ngs_pipeline_lib.tools.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external

from src.allele_issues import AlleleIssues
from src.inputs import MLSTInputs, SchemeMappingKB
from src.organism_comparison import ComparisonResult, compare_organism
from src.outputs import MLSTOutputs

# defined in the MLST image
MLST_SCHEME_SPECIES_MAP_PATH = getenv(
    "MLST_SCHEME_SPECIES_MAP_PATH", "/opt/conda/db/scheme_species_map.tab"
)
#MLST_SCHEME_SPECIES_MAP_PATH = "/scicomp/home-pure/lsz0/PN2.0/NewStructure/MLST/rdsi-ngs-pipeline-process-mlst/scheme_species_map.tab" ## no container

MINIMUM_COVERAGE_FOR_PARTIAL_ALLELES = 70


class MLST(Algorithm[MLSTInputs, MLSTOutputs]):
    outputs_class = MLSTOutputs
    working_dir = Path("mlst").resolve()

    def execute_stub(self):
        """
        Nothing to stub
        """
        pass

    def execute_implementation(self):
        """
        inputs:
        - assembled genome (preferably corrected assembly)
        outputs:
        - outputs.json: ST profile of the genome (if one can be found in the MLST schemes)
        - novel_fasta.tar.gz: an archive containing the novel_{scheme}.fasta in case there are any novel alleles not previously described
        """

        schemes = get_schemes_from_kb(
            self.inputs.scheme_mapping_kb, self.inputs.organism
        )

        self.result["mlst_results"] = []
        # Perfom the autodetection
        mlst_profile_json, novel_fasta = call_mlst(
            assembly=self.inputs.assembly,
            workdir=self.working_dir,
            logger=self.logger,
            n_threads=self.inputs.n_threads,
        )
        auto_detect_result = parse_mlst(mlst_profile_json)
        auto_detect_result["auto_detected"] = True
        if novel_fasta.exists():
            auto_detect_result["novel"] = str(novel_fasta) ## error fix: this is PosixPath which causes json error downstream, change to str
            if self.outputs.novel_fastas.content is not None: ## error fix: added check, cannot append before creating list
                self.outputs.novel_fastas.content.append(novel_fasta)
            else:
                self.outputs.novel_fastas.content = [novel_fasta]
        else:
            auto_detect_result["novel"] = None
        auto_detected_scheme = auto_detect_result["scheme"]
        self.result["mlst_results"].append(auto_detect_result)

        # For all the KB-defined schemes minus the one from auto-detection if it is in the list
        kb_schemes_result = []
        for scheme in schemes - set([auto_detected_scheme]):
            mlst_profile_json, novel_fasta = call_mlst(
                assembly=self.inputs.assembly,
                scheme=scheme,
                workdir=self.working_dir,
                logger=self.logger,
                n_threads=self.inputs.n_threads,
            )

            mlst_result = parse_mlst(mlst_profile_json)
            mlst_result["auto_detected"] = False
            if novel_fasta.exists():
                mlst_result["novel"] = str(novel_fasta) ## error fix: this is a PosixPath obj which causes json error downstream, change to str
                if self.outputs.novel_fastas.content is not None: ## error fix: added check, cannot append before creating list
                    self.outputs.novel_fastas.content.append(novel_fasta)
                else:
                    self.outputs.novel_fastas.content = [novel_fasta]
                
            else:
                mlst_result["novel"] = None
            kb_schemes_result.append(mlst_result)

        self.result["mlst_results"] += kb_schemes_result

        self.qc_report.add_metrics(
            compute_metrics(
                auto_detect_result=auto_detect_result,
                kb_schemes_result=kb_schemes_result,
                schemes=schemes,
                organism=self.inputs.organism,
            )
        )
        apply_quality_control(
            qc_report=self.qc_report,
            organism=self.inputs.organism,
            qc_dict=self.inputs.qc_kb.qc.get_dict(),
        )


def get_organisms_from_scheme_species_mapping(
    target_scheme: str,
) -> Iterable[OrganismInput]:
    with open(MLST_SCHEME_SPECIES_MAP_PATH, encoding="utf-8") as reader:
        for i, (scheme, genus, species) in enumerate(
            csv_reader(reader, delimiter="\t")
        ):
            if i == 0:
                continue
            if scheme.strip() == target_scheme.strip():
                yield OrganismInput(
                    genus=(genus.strip() or None), species=(species.strip() or None)
                )


def get_schemes_from_kb(
    mapping_kb: SchemeMappingKB, organism: OrganismInput
) -> set[str]:
    """
    Look inside the KB:
    if genus+species is found, return the schemes
    if not, look for genus and return all the schemes that have this genus,
    return an empty set if nothing matches the genus
    """
    if mapping_kb is None:
        return set()

    mapping_dict = mapping_kb.mapping.get_dict()

    try:
        species = organism.species if organism.species else "DEFAULT"
        return set(mapping_dict[organism.genus][species])
    except KeyError:
        result = set()
        for schemes in mapping_dict.get(organism.genus, {}).values():
            result |= set(schemes)
        return result


def call_mlst(
    assembly: Path,
    workdir: Path,
    logger: Logger,
    mincov: int | None = MINIMUM_COVERAGE_FOR_PARTIAL_ALLELES,
    scheme: str | None = None,
    n_threads: int | None = None,
) -> tuple[Path, Path]:
    """
    Calls [Torsten's MLST](https://github.com/tseemann/mlst) program to get the Sequence Type (ST) for a genome
    Args:
        assembly (Path): genome.fasta(.gz) that you want MLST to be run for
        logger (Logger): logger
        n_threads (int | None, optional): number of threads to use

    Returns:
        tuple[Path, Path]: Path to the JSON output, whether a novel.fasta was created
    """

    makedirs(workdir, exist_ok=True)
    mlst_profile_json = (workdir / "mlst_profile.json").resolve()

    # command_line = ["-n", "mlst", "mlst"]
    command_line = ["mlst"]
    command_line += [str(assembly.resolve())]
    if mincov is not None:
        command_line += ["--mincov", str(mincov)]
    if n_threads:
        command_line += ["--threads", str(n_threads)]
    if scheme is not None:
        novel_fasta = (workdir / f"novel_{scheme}.fasta").resolve()
        command_line += ["--scheme", scheme]
    else:
        novel_fasta = (workdir / "novel.fasta").resolve()
        command_line += ["--exclude", ""] ## do not exclude schemes by default during autodetect
    command_line += ["--nopath"]
    command_line += ["--json", str(mlst_profile_json)]
    command_line += ["--novel", str(novel_fasta)]

    _, stderr = run_external(
        command_line, logger, use_mamba_env=True, working_dir=workdir ## use_mamba_env=True for container, otherwise false
    )
    # The MLST program writes logs to STDERR, while STDOUT is reserved only for the output in TSV format
    # We will ignore STDOUT since all information in the TSV output is already present in the JSON file
    logger.info(stderr)
    return mlst_profile_json, novel_fasta


def parse_mlst(
    results: Path,
) -> dict[str, str | dict[str, str]]:
    """Simply parses over the JSON output to convert it into a dictionary

    Args:
        results (Path): path to the JSON output of the MLST program

    Returns:
        dict[str, str | dict[str, str]]: outputs the JSON content
    """
    with open(results, encoding="utf-8") as reader:
        output = json_load(reader)
    return output[0]


def get_best_organism_comparison_from_scheme_species_mapping(
    organism: OrganismInput, scheme: str
):
    """
    Checks all possible organisms a scheme can represent
    based on the scheme_species_map and then compares each
    of those organisms to the input organism, returning
    only the best possible comparison
    """
    observed_organisms = list(
        get_organisms_from_scheme_species_mapping(target_scheme=scheme)
    )
    return (
        min(
            compare_organism(
                observed_organism=observed_organism,
                process_organism=organism,
            )
            for observed_organism in observed_organisms
        )
        if observed_organisms
        else ComparisonResult.observed_organism_unknown
    )


def compute_metrics(
    auto_detect_result: dict[str, str | dict[str, str]],
    kb_schemes_result: list[dict[str, str | dict[str, str]]],
    schemes: list[str],
    organism: OrganismInput,
) -> dict[str, bool | str | int]:
    allele_issue_results: list[AlleleIssues] = []
    auto_detected_scheme = auto_detect_result["scheme"]

    auto_detected_organism_comparison = ComparisonResult.observed_organism_unknown

    if auto_detected_scheme in schemes:
        auto_detected_organism_comparison = ComparisonResult.ok
        kb_schemes_result.append(auto_detect_result)
    else:
        # if several organisms are linked to the scheme
        # we return the best result from the comparisons
        auto_detected_organism_comparison = (
            get_best_organism_comparison_from_scheme_species_mapping(
                organism, auto_detected_scheme
            )
        )
    auto_detected_allele_issues = AlleleIssues.from_scheme_result(
        auto_detect_result, auto_detected_organism_comparison
    )

    total_schemes = 0
    assigned_st = 0

    if not kb_schemes_result:
        allele_issues = auto_detected_allele_issues
        assigned_sequence_type_percent = float("nan")
        schemes_with_passing_alleles_percent = float("nan")
        schemes_with_allele_warnings_percent = float("nan")
        schemes_with_allele_errors_percent = float("nan")
    else:
        for scheme_result in kb_schemes_result:
            kb_result_organism_comparison = ComparisonResult.ok
            total_schemes += 1
            if scheme_result["sequence_type"] != "-":
                assigned_st += 1
            kb_result_allele_issues = AlleleIssues.from_scheme_result(
                scheme_result, organism_comparison=kb_result_organism_comparison
            )
            allele_issue_results.append(kb_result_allele_issues)
        allele_issues = (
            min(allele_issue_results)
            if allele_issue_results
            else auto_detected_allele_issues
        )
        # if several schemes were used here, we select the one with
        # the least amount of errors in it (not necessarily the best
        # result from comparisons)
        assigned_sequence_type_percent = assigned_st / (total_schemes or 1) * 100.0
        passing_schemes_issues = [
            scheme_allele_issues
            for scheme_allele_issues in allele_issue_results
            if not scheme_allele_issues.has_warning
            and not scheme_allele_issues.has_error
        ]
        schemes_with_passing_alleles_percent = (
            len(passing_schemes_issues) / (total_schemes or 1) * 100.0
        )
        schemes_with_allele_warnings_percent = (
            len(
                [
                    scheme_allele_issues
                    for scheme_allele_issues in allele_issue_results
                    if scheme_allele_issues.has_warning
                    and not scheme_allele_issues.has_error
                ]
            )
            / (total_schemes or 1)
            * 100.0
        )
        schemes_with_allele_errors_percent = (
            len(
                [
                    scheme_allele_issues
                    for scheme_allele_issues in allele_issue_results
                    if scheme_allele_issues.has_error
                ]
            )
            / (total_schemes or 1)
            * 100.0
        )
    return {
        "autoDetectFailed": auto_detected_allele_issues.has_warning,
        "autoDetectedScheme": auto_detected_allele_issues.scheme,
        "autoDetectedOrganismComparison": str(auto_detected_organism_comparison),
        "bestScheme": allele_issues.scheme,
        "organismComparison": str(allele_issues.organism_comparison),
        "assignedSequenceTypePercent": assigned_sequence_type_percent,
        "schemesWithoutIssuesPercent": schemes_with_passing_alleles_percent,
        "schemesWithWarningsPercent": schemes_with_allele_warnings_percent,
        "schemesWithErrorsPercent": schemes_with_allele_errors_percent,
        "unknownAlleles": len(allele_issues.unknown_alleles),
        "missingAlleles": len(allele_issues.missing_alleles),
        "multipleAlleles": len(allele_issues.multiple_alleles),
        "totalAlleleIssues": allele_issues.total_issues,
    }


def apply_quality_control(
    qc_report: QCReport,
    organism: OrganismInput,
    qc_dict: dict,
):
    quality_control = QualityControl(
        qc_dict=qc_dict, report=qc_report, organism=organism
    )

    quality_control.apply(
        section_name="mlstTaxonomy",
        observations=qc_report.metrics,
    )
