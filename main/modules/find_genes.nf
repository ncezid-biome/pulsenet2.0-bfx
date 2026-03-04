process FIND_GENES {

    publishDir "$params.publish_dir/${sample}/FIND_GENES", pattern: "find_genes.tsv", mode: "copy"

    input:
        tuple val(sample), path(assembly), val(organism), path(masterOut)

    output:
        tuple val(sample), path("outputs.json"), path("AMRFinderResultsRaw.json"), emit: master
        tuple val(sample), path("find_genes.tsv"), emit: result
        tuple val(sample), path("logs/messages.log")

    script:
    """
    species="\$(grep '"SPECIES": ' $masterOut | awk -F '"' '{print \$4}')"

    ngs-run \
    --sample-id $sample \
    --publish-dir . \
    --assembly $assembly \
    --organism.genus $organism \
    --organism.species \$species \

    """

    stub:
    """
    # ngs-run \
    # --sample-id $sample \
    # --publish-dir . \
    # --assembly $assembly \
    # --organism.genus $organism \
    # --stub

    touch outputs.json
    touch AMRFinderResultsRaw.json
    touch find_genes.tsv
    mkdir logs
    touch logs/messages.log
    """
}  

