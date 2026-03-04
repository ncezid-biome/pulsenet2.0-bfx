process ALLELE_NAMING {

    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"

    input:
        tuple val(sample), val(org_genus), path(allele_calls_xml), path(standard_calls), path(tempMasterOut, stageAs: "tempPipelineProcessOutputs.json")
        path(allele_cache_kb)
        val(allele_cache_ref)
        val(nomenclature_settings)
        path(masterOutScript)

    output:
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), path("PipelineProcessOutputs.json"), emit: masterOut
        tuple val(sample), path("logs/messages.log")

    script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --allele-calls-xml $allele_calls_xml \
    --allele-calls-profile $standard_calls \
    --allele-cache-kb.path $allele_cache_kb \
    --allele-cache-kb.cache $allele_cache_ref/${org_genus}/ \
    --nomenclature-settings $nomenclature_settings/${org_genus}/nomenclature_setting.json

    python3 $masterOutScript --process AlleleNaming --tempMaster $tempMasterOut --parseFile outputs.json

    """

    stub:
    """
    #ngs-run \
    #--sample-id $sample \
    #--publish-dir . \
    #--allele-calls-xml $allele_calls_xml \
    #--allele-calls-profile $standard_calls \
    #--allele-cache-kb.path $allele_cache_kb \
    #--allele-cache-kb.cache $allele_cache_ref \
    #--nomenclature-settings $nomenclature_settings \
    #--stub

    touch outputs.json
    mkdir logs
    touch logs/messages.log
    touch PipelineProcessOutputs.json
    """
}