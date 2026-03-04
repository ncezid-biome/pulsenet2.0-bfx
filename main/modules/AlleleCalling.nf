process ALLELE_CALLING {

    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"

    input:
        tuple val(sample), path(assembly), val(org_genus), path(tempMasterOut, stageAs: "tempPipelineProcessOutputs.json"), val(similarity)
        path(blast_kb)
        val(blastdb)
        val(loci)
        path(qc_kb)
        path(masterOutScript)

    output:
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), path("stats_calls.json.gz"), emit: stats
        tuple val(sample), path("allele_calls.bam"), path("PipelineProcessOutputs.json"), env(QC), emit: allele_calls_bam
        tuple val(sample), path("allele_calls.xml.gz"), emit: allele_calls_xml
        tuple val(sample), path("allele_calls.json.gz"), emit: allele_calls_json
        tuple val(sample), path("calls_standard.json.gz"), emit: standard_calls
        tuple val(sample), path("calls_core_standard.csv.gz"), path("calls_core_pcr.csv.gz"), emit: csv_core
        tuple val(sample), path("calls_accessory_standard.csv.gz"), path("calls_accessory_pcr.csv.gz"), emit: csv_accessory
        tuple val(sample), path("logs/messages.log")

    script:
    """

    species="\$(grep '"SPECIES": ' $tempMasterOut | awk -F '"' '{print \$4}')"
    lociList="${org_genus}_loci.csv"
    if [[ $org_genus == "VIBRIO" ]]; then
        if [[ "\$species" == "parahaemolyticus" || "\$species" == "cholerae" ]]; then
            lociList="${org_genus}_\${species}_loci.csv"
        fi
    fi

    ngs-run AlleleCalling \
    --sample-id $sample \
    --publish-dir . \
    --n-threads 1 \
    --assembly $assembly \
    --blast-kb.similarity $similarity \
    --blast-kb.path $blast_kb \
    --blast-kb.db $blastdb/$org_genus \
    --subspecies $tempMasterOut \
    --blast-kb.loci $loci/\$lociList \
    --qc-kb.path $qc_kb \
    --organism.genus $org_genus \
    --organism.species \$species \

    QC="true"
    qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
    if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
        QC="true"
    else
        if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
            QC="false"
        fi
    fi

    python3 $masterOutScript --process AlleleCalling --tempMaster $tempMasterOut --parseFile outputs.json

    """

    stub:
    """
    #ngs-run AlleleCalling \
    #--sample-id $sample \
    #--publish-dir . \
    #--n-threads 1 \
    #--assembly $assembly \
    #--blast-kb.similarity $similarity \
    #--blast-kb.path $blast_kb \
    #--blast-kb.db $blastdb/$org_genus \
    #--blast-kb.loci $loci/${org_genus}_loci.csv \
    #--qc-kb.path $qc_kb \
    #--organism.genus $org_genus \
    #--stub

    
    QC="true"
    touch outputs.json
    touch stats_calls.json.gz
    touch allele_calls.bam
    touch allele_calls.xml.gz
    touch allele_calls.json.gz
    touch calls_standard.json.gz
    touch calls_core_standard.csv.gz
    touch calls_accessory_standard.csv.gz
    touch calls_core_pcr.csv.gz
    touch calls_accessory_pcr.csv.gz
    mkdir logs
    touch logs/messages.log
    touch PipelineProcessOutputs.json

    """

}

process ALLELE_FILTERING {

    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"

    input:
        tuple val(sample), path(assembly), val(org_genus), path(cram_alignment), path(allele_calls_bam, stageAs: "input_allele_calls.bam"), path(tempMasterOut, stageAs: "tempMasterOut.json")
        path(filtering_kb)
        path(qc_kb)
        path(masterOutScript)

    output:
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), path("stats_calls.json.gz"), emit: stats
        tuple val(sample), path("allele_calls.bam"), path("allele_calls.json.gz"), emit: allele_calls
        tuple val(sample), path("calls_standard.json.gz"), path("PipelineProcessOutputs.json"), env(QC), emit: standard_calls
        tuple val(sample), path("PipelineProcessOutputs.json"), emit: master
        tuple val(sample), path("calls_core_standard.csv.gz"), path("calls_core_pcr.csv.gz"), emit: csv_core
        tuple val(sample), path("calls_accessory_standard.csv.gz"), path("calls_accessory_pcr.csv.gz"), emit: csv_accessory
        tuple val(sample), path("logs/messages.log")

    script:
    """
    species="\$(grep '"SPECIES": ' $tempMasterOut | awk -F '"' '{print \$4}')"
    
    ngs-run AlleleFiltering \
    --sample-id $sample \
    --publish-dir . \
    --n-threads 1 \
    --assembly $assembly \
    --alignment $cram_alignment \
    --calls-bam $allele_calls_bam \
    --filtering-kb.path $filtering_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $org_genus \
    --organism.species \$species \
    --subspecies $tempMasterOut \


    QC="true"
    qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
    if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
      QC="true"
    else
        if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
            QC="false"
        fi
    fi
    
    python3 $masterOutScript --process AlleleFilter --tempMaster $tempMasterOut --parseFile outputs.json
    """

    stub:
    """
    #ngs-run AlleleFiltering \
    #--sample-id $sample \
    #--publish-dir . \
    #--n-threads 1 \
    #--assembly $assembly \
    #--alignment $cram_alignment \
    #--calls-bam $allele_calls_bam \
    #--filtering-kb.path $filtering_kb \
    #--qc-kb.path $qc_kb \
    #--organism.genus $org_genus \
    #--stub

    QC="true"
    touch outputs.json
    touch stats_calls.json.gz
    touch allele_calls.bam
    touch allele_calls.json.gz
    touch calls_standard.json.gz
    touch calls_core_standard.csv.gz
    touch calls_accessory_standard.csv.gz
    touch calls_core_pcr.csv.gz
    touch calls_accessory_pcr.csv.gz
    mkdir logs
    touch logs/messages.log
    touch PipelineProcessOutputs.json
    """
}

