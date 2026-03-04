process MLST {

    input:
        tuple val(sample), path(assembly), val(org_genus)
        path(qc_kb)
        path(scheme_mapping_kb)
        path(mlstCC_script)
        path(cc_mapping)

    output:
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), path("novel_fastas.tar.gz"), emit: novel
        tuple val(sample), path("logs/messages.log")

    script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --n-threads 1 \
    --qc-kb.path $qc_kb \
    --scheme-mapping-kb.path $scheme_mapping_kb \
    --organism.genus $org_genus \
    --assembly $assembly \

    python3 $mlstCC_script --mlstOutput outputs.json --mappingFile $cc_mapping
    """

    stub:
    """
    #ngs-run \
    #--sample-id $sample \
    #--publish-dir . \
    #--n-threads 1 \
    #--qc-kb.path $qc_kb \
    #--scheme-mapping-kb.path $scheme_mapping_kb \
    #--organism.genus $org_genus \
    #--assembly $assembly \
    #--stub

    touch outputs.json
    touch novel_fastas.tar.gz
    mkdir logs
    touch logs/messages.log
    """
}