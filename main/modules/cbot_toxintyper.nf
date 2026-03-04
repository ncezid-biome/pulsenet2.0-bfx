process CBOT_TOXINTYPER {

    input:
        tuple val(sample), path(assembly), val(org_genus)
        path(cbot_toxintyping_db)

    output:
        tuple val(sample), path("outputs.json"), path("cbot_toxintyper.json"), emit: master
        tuple val(sample), path("logs/messages.log")

    script:
    """
    ngs-run CbotToxinTyping \
    --sample-id $sample \
    --publish-dir . \
    --assembly $assembly \
    --organism.genus $org_genus \
    --reference $cbot_toxintyping_db \
    """

    stub:
    """
    ngs-run CbotToxinTyping \
    --sample-id $sample \
    --publish-dir . \
    --assembly $assembly \
    --organism.genus $org_genus \
    --reference $cbot_toxintyping_db \
    --stub
    """
}
