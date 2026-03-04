FAKE_FASTP_REPORT = {
    "summary": {
        "fastp_version": "0.23.2",
        "sequencing": "paired end (301 cycles + 301 cycles)",
        "before_filtering": {
            "total_reads": 1201746,
            "total_bases": 320653185,
            "q20_bases": 262615025,
            "q30_bases": 226942024,
            "q20_rate": 0.819,
            "q30_rate": 0.707749,
            "read1_mean_length": 260,
            "read2_mean_length": 272,
            "gc_content": 0.392317,
        },
        "after_filtering": {
            "total_reads": 1097212,
            "total_bases": 284082366,
            "q20_bases": 242797536,
            "q30_bases": 213589710,
            "q20_rate": 0.854673,
            "q30_rate": 0.751858,
            "read1_mean_length": 257,
            "read2_mean_length": 260,
            "gc_content": 0.387109,
        },
    },
    "filtering_result": {
        "passed_filter_reads": 1097212,
        "low_quality_reads": 103126,
        "too_many_N_reads": 1408,
        "too_short_reads": 0,
        "too_long_reads": 0,
    },
    "duplication": {"rate": 0.000436032},
    "command": " ".join(
        [
            "fastp --detect_adapter_for_pe",
            "--in1 downsampled_read_1.fastq.gz",
            "--in2 downsampled_read_2.fastq.gz",
            "--out1 cleaned_read_1.fastq",
            "--out2 cleaned_read_2.fastq",
            "-j fastp_report.json ",
        ]
    ),
}
