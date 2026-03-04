process MIDAS {

  publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"
  publishDir "$params.publish_dir/${sample}/", pattern: "ADX_analysis_data.json", mode : "copy"

  input:
    tuple val(sample), path(read1), path(read2), path(tempMasterOut, stageAs: "tempMasterOut.json"), val(org_genus)
    val(midas)
    path(qc_kb)
    val(process)
    path(masterOutScript)

  output:
    tuple val(sample), path("outputs.json"), emit: outputs
    tuple val(sample), path("midas_result_out.json"), emit: midas_result
    tuple val(sample), path("species/species_profile.txt"), emit: species_profile
    tuple val(sample), stdout, path("PipelineProcessOutputs.json"), env(QC), emit: genus
    tuple val(sample), path("logs/messages.log"), path("ADX_analysis_data.json")

  script:
    """
    echo "[]" > ADX_analysis_data.json

    export MIDAS_DB=${midas}
    ngs-run \
    --publish-dir . \
    --sample-id $sample \
    --read1 $read1 \
    --read2 $read2 \
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

    if [ $process == "GENUS_IDENTIFICATION" ]; then
      python3 $masterOutScript --process GenusIdentify --parseFile outputs.json
    else
      python3 $masterOutScript --process Midas --tempMaster $tempMasterOut --parseFile outputs.json --parseAddition midas_result_out.json
    fi
    """

 
  stub:
    """
    touch outputs.json
    touch midas_result_out.json
    mkdir species
    touch species/species_profile.txt
    touch PipelineProcessOutputs.json
    touch temp_PipelineProcessOutputs.json
    mkdir logs
    touch logs/messages.log
    QC="true"

    genus[0]="Salmonella"
    genus[1]="Escherichia"
    genus[2]="Listeria"
    #genus[3]="Vibrio"
    #genus[4]="Campylobacter"
    #genus[5]="Cronobacter"
    #genus[6]="Yersinia"
    #genus[7]="Clostridium"

    sample=\$(cksum <(echo $sample) | cut -d ' ' -f1)
    i=\$(( \$sample % \${#genus[@]} ))
    printf "\${genus[\$i]}"
    """
}
