process GENERATE_ASSEMBLY {

  publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"

  input:
    tuple val(sample), path(read1), path(read2), val(org_genus), path(tempMasterOut, stageAs: "tempMasterOut.json")
    path(qc_kb)
    path(masterOutScript)

  output:
    tuple val(sample), path("outputs.json"), emit: outputs
    tuple val(sample), path("assembly.fasta.gz"), path("PipelineProcessOutputs.json"), env(QC), emit: assembly
    tuple val(sample), path("logs/messages.log")

  script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --read1 $read1 \
    --read2 $read2 \
    --qc-kb.path $qc_kb \
    --organism.genus $org_genus \
    --n-threads ${task.cpus}

    QC="true"
    qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
    if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
      QC="true"
    else
      if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
          QC="false"
      fi
    fi
    
    python3 $masterOutScript --process GenerateAssembly --tempMaster $tempMasterOut --parseFile outputs.json

    """

  stub:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --read1 $read1 \
    --read2 $read2 \
    --qc-kb.path $qc_kb \
    --organism.genus $org_genus \
    --n-threads ${task.cpus} \
    --stub

    QC="true"
    touch PipelineProcessOutputs.json
    """
}