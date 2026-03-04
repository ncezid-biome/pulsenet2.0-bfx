process STXCONDENSER {
   
    input:
        tuple val(sample), path(stx_type_result), path(isprc_result)
    output:
        tuple val(sample), path("stx_condenser_expr.json"), path("stx_condenser_flds.json"), emit: stxcondenser_results
        tuple val(sample), path("outputs.json"), path("master_input.json"), emit: master
        tuple val(sample), path("stx_condenser_expr.json"), path("stx_condenser_flds.json"),
        path ("*"), emit: allOut
        tuple val(sample), path("logs/messages.log")

    script:
        """
        ngs-run \
        --sample-id $sample \
        --publish-dir ./ \
        --stx-results $stx_type_result \
        --ispcr-results $isprc_result

        echo "[\$(cat stx_condenser_expr.json),\$(cat stx_condenser_flds.json)]" > master_input.json
        """

    stub:
        """
        touch stx_condenser_expr.json
        touch stx_condenser_flds.json
        touch outputs.json
        mkdir logs
        touch logs/messages.log
        """
}
