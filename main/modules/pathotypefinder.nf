process PATHOTYPE_FINDER {

    input:
        tuple val(sample), path(assembly), val(org_genus)
        path(pathotype_db)

    output:
        tuple val(sample), path("outputs.json"), path("pathotypefinder_genotypes.json"), emit: master
        tuple val(sample), path("blastout.txt")
        tuple val(sample), path("logs/messages.log")

    script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --assembly $assembly \
    --organism.genus $org_genus \
    --reference $pathotype_db \
    """

    stub:
    """
    # ngs-run \
    # --sample-id $sample \
    # --publish-dir . \
    # --assembly $assembly \
    # --organism.genus $org_genus \
    # --reference $pathotype_db \
    # --stub
    touch pathotypefinder_genotypes.json
    touch blastout.txt
    touch outputs.json
    mkdir logs
    touch logs/messages.log
    """
}