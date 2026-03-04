import os
from abc import ABCMeta, abstractmethod
from collections import Counter
from dataclasses import dataclass
from functools import partial
from logging import Logger
from pathlib import Path
from shutil import rmtree
from typing import Iterable

from Bio.SeqIO import parse, to_dict, write  # type: ignore
from Bio.SeqRecord import SeqRecord
from ngs_pipeline_lib.tools.parallel import run_parallel
from ngs_pipeline_lib.tools.runextern import run_external

##CONDA_BIN_DIR = os.getenv("CONDA_BIN_DIR", "/opt/conda/condabin") #Container


@dataclass
class CleanedAssemblyInfo:
    initial_assembly_length: int
    initial_average_depth: float
    initial_contigs_number: int
    initial_gc_content: float
    initial_general_horizontal_coverage: int
    initial_n50: int

    cleaned_assembly_length: int
    cleaned_average_depth: float
    cleaned_contigs_number: int
    cleaned_gc_content: float
    cleaned_general_horizontal_coverage: int
    cleaned_n50: int


class ContigInfo(metaclass=ABCMeta):
    def __init__(
        self, contigs_dict: dict[str, SeqRecord], contig_list: None | list[str] = None
    ) -> None:
        if contig_list is None:
            self.contig_list = contigs_dict.keys()
        else:
            self.contig_list = contig_list
        self.number_contigs = len(self.contig_list)

    @abstractmethod
    def get_message(self, prefix) -> str:
        ...

    @property
    @abstractmethod
    def metric(self):
        ...


class LengthContigInfo(ContigInfo):
    def __init__(
        self, contigs_dict: dict[str, SeqRecord], contig_list: None | list[str] = None
    ) -> None:
        super().__init__(contigs_dict, contig_list)
        self.length_contigs = list(
            map(lambda x: len(contigs_dict[x].seq), self.contig_list)
        )
        self.total_bp = sum(self.length_contigs)
        if self.number_contigs > 0:
            self.average_bp = round(self.total_bp / self.number_contigs, 2)
        else:
            self.average_bp = 0

        ##
        self.n50 = 0
        n = self.total_bp / 2
        cumulated_length = 0
        for length in self.length_contigs:
            cumulated_length += length
            if cumulated_length >= n:
                self.n50 = length
                break
        ##

    def get_message(self, prefix: str):
        min_bp = min(self.length_contigs, default=0)
        max_bp = max(self.length_contigs, default=0)
        return (
            f"{prefix} {self.number_contigs=}; {self.total_bp=}; "
            f"{self.average_bp=}; {min_bp=}; {max_bp=}"
        )

    @property
    def metric(self):
        return self.total_bp


class GCContigInfo(ContigInfo):
    def __init__(
        self,
        contigs_dict: dict[str, SeqRecord],
        contig_list: None | Iterable[str] = None,
    ) -> None:
        super().__init__(contigs_dict, contig_list)
        self.length_contigs = list(
            map(lambda x: len(contigs_dict[x].seq), self.contig_list)
        )
        self.gc_contigs = list(
            map(
                lambda x: sum(
                    v
                    for k, v in Counter(contigs_dict[x].seq.upper()).items()
                    if k in ("G", "C")
                ),
                self.contig_list,
            )
        )
        self.gc_content = (
            round(sum(self.gc_contigs) / sum(self.length_contigs) * 100, 2)
            if self.number_contigs > 0
            else 0
        )
        self.gc_content_contigs = [
            round(v / self.length_contigs[x] * 100, 2)
            for x, v in enumerate(self.gc_contigs)
        ]

    def get_message(self, prefix: str):
        min_gc = min(self.gc_content_contigs, default=0)
        max_gc = max(self.gc_content_contigs, default=0)
        return (
            f"{prefix} {self.number_contigs=}; {self.gc_content=}; {min_gc=}; {max_gc=}"
        )

    @property
    def metric(self):
        return self.gc_content


class DepthContigInfo(ContigInfo):
    def __init__(
        self,
        depth_dict: dict[str, dict[str, float]],
        contig_list: None | list[str] = None,
    ) -> None:
        super().__init__(depth_dict, contig_list)
        self.depth_contigs = list(
            map(
                lambda x: round(depth_dict[x]["vertical"], 2),
                self.contig_list,
            )
        )
        try:
            self.average_depth = round(sum(self.depth_contigs) / self.number_contigs, 2)
        except ZeroDivisionError:
            self.average_depth = 0

    def get_message(self, prefix: str):
        min_avg_depth = min(self.depth_contigs, default=0)
        max_avg_depth = max(self.depth_contigs, default=0)
        return (
            f"{prefix} {self.number_contigs=}; {self.average_depth=}"
            f"; {min_avg_depth=} ; {max_avg_depth=}"
        )

    @property
    def metric(self):
        return None


@dataclass
class CleaningSettings:
    roles: list[str] = ("gc", "length", "depth")
    min_base_depth: int = 10
    min_depth: int = 15
    min_depth_percent: int = 25
    min_gc_content: int = 5
    min_contig_length: int = 100 ##added


def get_organism_cleaning_settings(
    organism_to_settings: dict[str, dict], genus, species
) -> CleaningSettings:
    settings = CleaningSettings()
    settings_dict = organism_to_settings.get(
        f"{genus} {species}",
        organism_to_settings.get(genus, organism_to_settings.get("DEFAULT", {})),
    )
    for setting, value in settings_dict.items():
        settings.__setattr__(setting, value)
    return settings


def log_contig_infos(
    logger: Logger,
    contigs_dict: dict[str, SeqRecord],
    depth_dict: dict[str, dict[str:float]],
    contig_list: list[str],
    prefix: str,
):
    length_contig_info = LengthContigInfo(contigs_dict, contig_list)
    logger.info(length_contig_info.get_message(prefix))

    gc_contig_info = GCContigInfo(contigs_dict, contig_list)
    logger.info(gc_contig_info.get_message(prefix))

    depth_contig_info = DepthContigInfo(depth_dict, contig_list)
    logger.info(depth_contig_info.get_message(prefix))


def process_samtools_stats_stdout(stdout: str, *metrics: str) -> tuple[float, ...]:
    results = {metric: 0 for metric in metrics}
    for line in stdout.split("\n"):
        stripped_line = line.strip()
        if stripped_line.startswith("SN"):
            split_line = stripped_line.split("\t")
            for metric in metrics:
                if split_line[1] == metric:
                    results[metric] = float(split_line[2])
    return tuple(results.values())


def compute_stats_for_contig(
    index: int,
    contig_name: str,
    contig: SeqRecord,
    fasta_to_use: Path,
    input_alignment: Path,
    align_index: Path,
    min_base_depth: int,
    logger: Logger,
):
    contig_len = len(contig.seq)
    target_regions_tsv = f"target-regions.{index}.tsv"
    with open(target_regions_tsv, "w", encoding="utf-8", newline="\n") as writer:
        print("\t".join(map(str, [contig_name, 1, contig_len])), file=writer)

    ##cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "--no-capture-output", "samtools", "stats"] #Container
    cmd = ["mamba", "run", "--no-capture-output", "samtools", "stats"] #Scicomp
    cmd += ["--ref-seq", str(fasta_to_use)]
    cmd += ["--target-regions", target_regions_tsv]
    cmd += ["--cov-threshold", str(min_base_depth - 1)]
    cmd += ["-@", "0"]
    cmd += ["-X", str(input_alignment)]
    cmd += [str(align_index)]

    stats_process_stdout, _ = run_external(cmd, logger)

    bases_mapped, horizontal_coverage = process_samtools_stats_stdout(
        stats_process_stdout,
        "bases mapped (cigar):",
        f"percentage of target genome with coverage > {min_base_depth - 1} (%):",
    )
    vertical_coverage = bases_mapped / contig_len
    return contig_name, contig_len, vertical_coverage, horizontal_coverage


def filter_contigs_on_depth(
    contigs_dict: dict[str, SeqRecord],
    min_contig_depth: float,
    min_base_depth: int,
    fasta_to_use: Path,
    input_alignment: Path,
    align_index: Path,
    workdir: Path,
    n_threads: int,
    logger: Logger,
) -> tuple[dict, dict]:
    depth_dict = {}
    output_dict = {"contigs": [], "length": []}
    depth_dir = workdir / "depth"
    depth_dir.mkdir(exist_ok=True)
    inputs = [
        {"index": index, "contig_name": contig_name, "contig": contig}
        for index, (contig_name, contig) in enumerate(contigs_dict.items())
    ]
    compute_stats_for_contig_with_files = partial(
        compute_stats_for_contig,
        min_base_depth=min_base_depth,
        fasta_to_use=fasta_to_use,
        input_alignment=input_alignment,
        align_index=align_index,
    )
    results = run_parallel(
        func=compute_stats_for_contig_with_files,
        inputs=inputs,
        n_threads=n_threads,
        logger=logger,
    )
    # to keep the same order as before (the multiprocessing does not keep it)
    results_as_dict = {result[0]: result[1:] for result in results}
    for contig_name in contigs_dict.keys():
        contig_len, vertical_coverage, horizontal_coverage = results_as_dict[
            contig_name
        ]
        depth_dict[contig_name] = {
            "vertical": vertical_coverage,
            "horizontal": horizontal_coverage,
        }

        if (vertical_coverage) < min_contig_depth:
            output_dict["contigs"].append(contig_name)
            output_dict["length"].append(contig_len)
    depth_contig_info = DepthContigInfo(depth_dict, None)
    logger.info(depth_contig_info.get_message("[CONTIG-DEPTH]"))

    if len(output_dict["contigs"]) > 0:
        log_contig_infos(
            logger,
            contigs_dict,
            depth_dict,
            output_dict["contigs"],
            "[FILTERED-CONTIG-DEPTH]",
        )
    rmtree(depth_dir)
    return depth_dict, output_dict


def filter_contigs_on_length(
    contigs_dict: dict[str, SeqRecord],
    depth_dict: dict[str, dict[str, float]],
    min_contig_length: int,
    logger: Logger,
) -> dict:
    output_dict = {"contigs": [], "length": []}
    for contig_name, contig_data in contigs_dict.items():
        contig_len = len(contig_data.seq)
        if contig_len < min_contig_length:
            output_dict["contigs"].append(contig_name)
            output_dict["length"].append(contig_len)

    if len(output_dict["contigs"]) > 0:
        log_contig_infos(
            logger,
            contigs_dict,
            depth_dict,
            output_dict["contigs"],
            "[FILTERED-CONTIG-LENGTH]",
        )
    return output_dict


def filter_contigs_on_gc(
    contigs_dict: dict[str, SeqRecord],
    depth_dict: dict[str, dict[str, float]],
    min_contig_gc,
    logger: Logger,
) -> dict:
    output_dict = {"contigs": [], "length": []}
    for contig_name, contig_data in contigs_dict.items():
        contig_len = len(contig_data.seq)
        base_counter = Counter(contig_data.seq.upper())
        gc_counts = sum(v for k, v in base_counter.items() if k in ("G", "C"))
        if min_contig_gc <= gc_counts / contig_len * 100 <= 100 - min_contig_gc:
            continue
        else:
            output_dict["contigs"].append(contig_name)
            output_dict["length"].append(contig_len)

    if len(output_dict["contigs"]) > 0:
        log_contig_infos(
            logger,
            contigs_dict,
            depth_dict,
            output_dict["contigs"],
            "[FILTERED-CONTIG-GC]",
        )
    return output_dict


def compute_stats_with_target_regions(
    contigs_dict: dict[str, SeqRecord],
    contig_list: Iterable[str],
    input_fasta: Path,
    input_alignment: Path,
    align_index: Path,
    workdir: Path,
    min_base_depth: int,
    n_threads: int,
    logger: Logger,
) -> str:
    # Create a 1-based index map of the coordinates from the entire contig
    # a list of genes coordinates, where a gene (or any genomic feat) starts and ends
    target_regions_tsv = workdir / "target-regions.assembly.tsv"
    with open(target_regions_tsv, "w", encoding="utf-8", newline="\n") as writer:
        for contig_name in contig_list:
            contig_len = len(contigs_dict[contig_name].seq)
            print("\t".join(map(str, [contig_name, 1, contig_len])), file=writer)

    ##cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "--no-capture-output", "samtools", "stats"] #Container
    cmd = ["mamba", "run", "--no-capture-output", "samtools", "stats"] #Scicomp
    cmd += ["--ref-seq", str(input_fasta)]
    cmd += ["-@", str(n_threads - 1)]
    cmd += ["--target-regions", str(target_regions_tsv)]
    cmd += ["--cov-threshold", str(min_base_depth - 1)]
    cmd += ["-X", str(input_alignment)]
    cmd += [str(align_index)]

    stats_process_stdout, _ = run_external(cmd, logger)
    target_regions_tsv.unlink()

    return stats_process_stdout


def create_bam_file(
    contigs_dict: dict[str, SeqRecord],
    cleaned_contigs: Iterable[str],
    bam: Path,
    input_fasta: Path,
    input_alignment: Path,
    align_index: Path,
    workdir: Path,
    n_threads: int,
    logger: Logger,
):
    target_regions = workdir / "target-regions.cleaned_contigs.bed"
    with open(
        target_regions,
        "w",
        encoding="utf-8",
        newline="\n",
    ) as writer:
        for contig_name in cleaned_contigs:
            print(
                "\t".join(
                    map(str, [contig_name, 0, len(contigs_dict[contig_name].seq)])
                ),
                file=writer,
            )
    ##cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "--no-capture-output", "samtools", "view"] #Container
    cmd = ["mamba", "run", "--no-capture-output", "samtools", "view"] #Scicomp
    cmd += ["-b"]
    cmd += ["-h"]
    cmd += ["-o", str(bam)]
    cmd += ["--reference", str(input_fasta)]
    cmd += ["-M"]
    cmd += ["-L", str(target_regions)]
    cmd += ["-@", str(n_threads - 1)]
    cmd += [f"{input_alignment}##idx##{align_index}"]

    run_external(cmd, logger)
    target_regions.unlink()


def create_fasta_file(
    filepath: Path, contigs_dict: dict[str, SeqRecord], cleaned_contigs: Iterable[str]
):
    with open(filepath, "w", encoding="utf-8", newline="\n") as writer:
        _ = write(
            [contigs_dict[contig_name] for contig_name in cleaned_contigs],
            writer,
            "fasta",
        )


def create_depth_contigs_file(
    depth_contigs: Path,
    contigs_dict: dict[str, SeqRecord],
    cleaned_contigs: Iterable[str],
    depth_dict: dict[str, dict[str, float]],
):
    with open(depth_contigs, "w", encoding="utf-8", newline="\n") as writer:
        print(
            "\t".join(
                [
                    "#contig",
                    "len",
                    "depth_x",
                    "horizontal_coverage_percent",
                    "filtered_status",
                ]
            ),
            file=writer,
        )
        for contig_name, x in depth_dict.items():
            print(
                "\t".join(
                    map(
                        str,
                        [
                            contig_name,
                            len(contigs_dict[contig_name].seq),
                            round(x["vertical"], 2),
                            x["horizontal"],
                            "passed"
                            if contig_name in set(cleaned_contigs)
                            else "filtered",
                        ],
                    )
                ),
                file=writer,
            )


def assembly_cleaning(
    input_fasta: Path,
    output_fasta: Path,
    input_alignment: Path,
    output_bam: Path,
    output_depth_stats: Path,
    cleaning_settings: CleaningSettings,
    n_threads: int,
    logger: Logger,
) -> CleanedAssemblyInfo:
    """
    Cleans and provides information about the assembly pre and post cleaning
    """
    workdir = Path("assembly_cleaning_workdir")
    workdir.mkdir(parents=True, exist_ok=True)
    with open(input_fasta, encoding="utf-8") as reader:
        contigs: dict[str, SeqRecord] = to_dict(parse(reader, "fasta"))

    initial_length_contig_info = LengthContigInfo(contigs, None)
    initial_assembly_length = initial_length_contig_info.metric
    logger.info(initial_length_contig_info.get_message("[CONTIG]"))

    gc_contig_info = GCContigInfo(contigs)
    initial_gc_content = gc_contig_info.metric
    logger.info(gc_contig_info.get_message("[CONTIG]"))

    # Create index for the alignment file (bam) -> align_index
    align_index = workdir / input_alignment.with_suffix(".index").name
    ##cmd = [f"{CONDA_BIN_DIR}/mamba", "run", "--no-capture-output", "samtools", "index"] #Container
    cmd = ["mamba", "run", "--no-capture-output", "samtools", "index"] #Scicomp
    cmd += ["-@", str(n_threads - 1)]
    cmd += [str(input_alignment)]
    cmd += [str(align_index)]
    run_external(cmd, logger)

    stats_process_stdout = compute_stats_with_target_regions(
        contigs,
        contigs.keys(),
        input_fasta,
        input_alignment,
        align_index,
        workdir,
        cleaning_settings.min_base_depth,
        n_threads,
        logger,
    )
    (
        bases_mapped,
        max_length,
        insert_length,
        initial_general_horizontal_coverage,
    ) = process_samtools_stats_stdout(
        stats_process_stdout,
        "bases mapped (cigar):",
        "maximum length:",
        "insert size average:",
        f"percentage of target genome with coverage > {cleaning_settings.min_base_depth - 1} (%):",
    )

    ##min_contig_length = int(round(max_length + ((insert_length - max_length) / 2), 0)) ## calculating min contig length to filter by
    min_contig_length = cleaning_settings.min_contig_length

    logger.info(
        (
            f"[MIN-CONTIG-LEN] {min_contig_length} ({max_length=} ; {insert_length=})"
            " (max_length + ((insert_length - max_length) / 2))"
        )
    )

    initial_average_depth = round(
        bases_mapped / sum(map(lambda seq_record: len(seq_record.seq), contigs.values())),
        2
        )
    logger.info(
        (
            f"[DEPTH] average {round(initial_average_depth, 2)} x, "
            f"{initial_general_horizontal_coverage}% covered at minimum of {cleaning_settings.min_base_depth}"
        )
    )

    min_depth = max(
        [
            cleaning_settings.min_depth,
            round(
                initial_average_depth * (cleaning_settings.min_depth_percent / 100), 2
            ),
        ]
    )
    # filering the contigs
    filtered_contigs = {}
    #   on depth
    depth_dict, depth_filtered_contigs = filter_contigs_on_depth(
        contigs,
        min_depth,
        cleaning_settings.min_base_depth,
        input_fasta,
        input_alignment,
        align_index,
        workdir,
        n_threads,
        logger,
    )
    filtered_contigs["depth"] = depth_filtered_contigs
    #   on length
    filtered_contigs["length"] = filter_contigs_on_length(
        contigs, depth_dict, min_contig_length, logger
    )
    #   on GC
    filtered_contigs["gc"] = filter_contigs_on_gc(
        contigs, depth_dict, cleaning_settings.min_gc_content, logger
    )

    # final cleaning
    filtered_contigs_final = set(
        contig_name
        for role in cleaning_settings.roles
        for contig_name in filtered_contigs[role]["contigs"]
    )

    log_contig_infos(
        logger, contigs, depth_dict, filtered_contigs_final, "[FILTERED-CONTIG]"
    )

    # subtract the filtered contigs but keep the order
    cleaned_contigs = [
        contig_name
        for contig_name in contigs.keys()
        if contig_name not in filtered_contigs_final
    ]

    cleaned_length_contig_info = LengthContigInfo(contigs, cleaned_contigs)
    cleaned_assembly_length = cleaned_length_contig_info.metric
    logger.info(cleaned_length_contig_info.get_message("[CLEANED-CONTIG]"))

    gc_contig_info = GCContigInfo(contigs, cleaned_contigs)
    cleaned_gc_content = gc_contig_info.metric
    logger.info(gc_contig_info.get_message("[CLEANED-CONTIG]"))

    depth_contig_info = DepthContigInfo(depth_dict, cleaned_contigs)
    logger.info(depth_contig_info.get_message("[CLEANED-CONTIG]"))

    stats_process_stdout = compute_stats_with_target_regions(
        contigs,
        cleaned_contigs,
        input_fasta,
        input_alignment,
        align_index,
        workdir,
        cleaning_settings.min_base_depth,
        n_threads,
        logger,
    )
    bases_mapped, cleaned_general_horizontal_coverage = process_samtools_stats_stdout(
        stats_process_stdout,
        "bases mapped (cigar):",
        f"percentage of target genome with coverage > {cleaning_settings.min_base_depth - 1} (%):",
    )

    if bases_mapped > 0:
        cleaned_average_depth = round(
            bases_mapped / sum(map(lambda x: len(contigs[x].seq), cleaned_contigs)),
            2,
        )
    else:
        cleaned_average_depth = 0

    create_fasta_file(output_fasta, contigs, cleaned_contigs)
    create_bam_file(
        contigs,
        cleaned_contigs,
        output_bam,
        input_fasta,
        input_alignment,
        align_index,
        workdir,
        n_threads,
        logger,
    )
    create_depth_contigs_file(output_depth_stats, contigs, cleaned_contigs, depth_dict)

    input_fasta.unlink(missing_ok=True)
    input_fasta.with_suffix(input_fasta.suffix + ".fai").unlink()
    align_index.unlink(missing_ok=True)

    generated_files = Path().glob("assembly.bam*")
    for generated_file in generated_files:
        generated_file.unlink()
    # return the data that will be used for the report.json
    return CleanedAssemblyInfo(
        initial_assembly_length,
        initial_average_depth,
        len(contigs),
        initial_gc_content,
        initial_general_horizontal_coverage,
        initial_length_contig_info.n50,
        cleaned_assembly_length,
        cleaned_average_depth,
        len(cleaned_contigs),
        cleaned_gc_content,
        cleaned_general_horizontal_coverage,
        cleaned_length_contig_info.n50,
    )
