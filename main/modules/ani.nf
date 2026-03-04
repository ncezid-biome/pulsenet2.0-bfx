process ANI{

    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"
    publishDir "$params.publish_dir/${sample}/", pattern: "ADX_analysis_data.json", mode : "copy"

    input:
        path(qc)
        tuple val(sample), path(assembly), path(tempMasterOut, stageAs: "tempMasterOut.json")
        val ref
        val emit
        path(masterOutScript)

    output:
        tuple val(sample), path("out.tsv"), emit: all
        tuple val(sample), path("best_hit.json"), emit: best
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), stdout, path("PipelineProcessOutputs.json"), env(QC), emit: organism
        tuple val(sample), path("logs/messages.log")
        tuple val(sample), path("*.json")

    script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir ./ \
    --assembly $assembly \
    --reference $ref \
    --qc $qc/qc.json \
    --qc-section "ANI" \
    --n-threads $task.cpus 
    
    python3 $masterOutScript --process ANI --tempMaster $tempMasterOut --parseFile outputs.json --parseAddition best_hit.json

    QC="true"
    qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
    if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
            QC="true"
    else
        if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
            QC="false"
        fi
    fi

    #Check if PulseNet Genus and Species
    pn_genus="Salmonella Escherichia Campylobacter Listeria Cronobacter Vibrio Clostridium Yersinia"
    genus_to_check="\$(grep '"GENUS": ' PipelineProcessOutputs.json | awk -F '"' '{print \$4}')"

    if ! \$(echo "\$pn_genus" | grep -qw "\$genus_to_check"); then
        echo "[]" > ADX_analysis_data.json
        QC="false"
    else
        if [[ $params.non_surv_orgs == "false" ]]; then
            pn_species="botulinum cholerae coli enterica enterocolitica fetus jejuni lari monocytogenes parahaemolyticus sakazakii upsaliensis vulnificus albertii condimenti dublinensis malonaticus muytjensii turicensis universalis"
            species_to_check="\$(grep '"SPECIES": ' PipelineProcessOutputs.json | awk -F '"' '{print \$4}')"

            if ! \$(echo "\$pn_species" | grep -qw "\$species_to_check"); then
                echo "[]" > ADX_analysis_data.json
                QC="false"
            fi
        fi
    fi
   
    """

    stub:
    """
    # ngs-run \
    # --sample-id $sample \
    # --publish-dir ./ \
    # --assembly $assembly \
    # --reference $ref \
    # --qc $qc \
    # --n-threads $task.cpus \
    # --stub

    touch out.tsv
    touch best_hit.json
    touch outputs.json
    touch PipelineProcessOutputs.json
    mkdir logs
    touch logs/messages.log

    QC="true"

    genus[0]="SALM"
    genus[1]="STEC"
    genus[2]="LISTERIA"
    
    #genus[3]="VIBRIO"
    #genus[4]="CAMPY"
    #genus[5]="CRONO"
    #genus[6]="CBOT"
    #genus[7]="YERSINIA"

    sample=\$(cksum <(echo $sample) | cut -d ' ' -f1)
    i=\$(( \$sample % \${#genus[@]} ))
    printf "\${genus[\$i]}"
    """
}

process ORG_ANI{

    input:
        path(qc)
        tuple val(sample), path(assembly), val(org)
        val(salm_ref)
        val emit

    output:
        tuple val(sample), path("out.tsv"), emit: all
        tuple val(sample), path("best_hit.json"), env(QC), emit: best
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), stdout, emit: organism
        tuple val(sample), path("logs/messages.log")

    script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir ./ \
    --assembly $assembly \
    --reference $salm_ref \
    --qc $qc/qc.json \
    --qc-section "ORG_ANI" \
    --n-threads $task.cpus 

    QC="true"
    qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
    if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
      QC="true"
    else
        if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
            QC="false"
        fi
    fi

    """

    stub:
    """
    # ngs-run \
    # --sample-id $sample \
    # --publish-dir ./ \
    # --assembly $assembly \
    # --reference $salm_ref \
    # --qc $qc \
    # --n-threads $task.cpus \
    # --stub

    touch out.tsv
    touch best_hit.json
    touch outputs.json
    mkdir logs
    touch logs/messages.log

    QC="true"
    printf "SALM"
    """


}