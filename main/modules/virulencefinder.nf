process VIRULENCE_FINDER {

    input:
        tuple val(sample), path(assembly), val(organism)
        val(virulencefinder_db)

    output:
        tuple val(sample), path("outputs.json"), path("virulence.json"), emit: master
        tuple val(sample), path("blastout.txt")
        tuple val(sample), path("logs/messages.log")

    script:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --assembly $assembly \
    --organism.genus $organism \
    --reference $virulencefinder_db \
    """

    stub:
    """
    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --assembly $assembly \
    --organism.genus $org_genus \
    --reference $virulencefinder_db \
    --stub
    """
}
