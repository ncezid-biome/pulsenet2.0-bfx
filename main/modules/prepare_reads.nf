process PREPARE_READS {

  publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"

  input:
    tuple val(sample), path(read1), path(read2), val(org_genus), path(tempMasterOut, stageAs: "tempMasterOut.json")
    path(genome_sizes)
    path(qc_kb)
    path(masterOutScript)

  output:
    tuple val(sample), path("outputs.json"), emit: outputs
    tuple val(sample), path("cleaned_read_1.fastq.gz"), path("cleaned_read_2.fastq.gz"), path("PipelineProcessOutputs.json"), env(QC), emit: reads
    tuple val(sample), path("fastp_report.json"), emit: fastp
    tuple val(sample), path("logs/messages.log")

  script:
    """
    ngs-run CleanupReads \
    --sample-id $sample\
    --publish-dir . \
    --read1 $read1 \
    --read2 $read2 \
    --genome-kb.path $genome_sizes \
    --qc-kb.path $qc_kb \
    --organism.genus $org_genus

    QC="true"
    qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
    if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
      QC="true"
    else
      if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
          QC="false"
      fi
    fi
    python3 $masterOutScript --process PrepareReads --tempMaster $tempMasterOut --parseFile outputs.json

    """

  stub:
    """
    ngs-run CleanupReads \
    --sample-id $sample\
    --publish-dir . \
    --read1 $read1 \
    --read2 $read2 \
    --genome-kb.path $genome_sizes \
    --qc-kb.path $qc_kb \
    --organism.genus $org_genus \
    --stub

    QC="true"
    touch PipelineProcessOutputs.json

    """
}
