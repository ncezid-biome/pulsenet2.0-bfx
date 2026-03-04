process PLASMIDFINDER {

  input:
       	tuple val(sample), path(assembly), val(organism)
        val(ref)
  output:
        tuple val(sample), path("outputs.json"), path("plasmidfinder.json"), emit: master
        path("*")
        tuple val(sample), path("logs/messages.log")

  script:
    """
    ngs-run \
    --publish-dir . \
    --sample-id $sample \
    --assembly $assembly \
    --organism.genus $organism \
    --reference $ref \ 
    """

  stub:
    """
    # ngs-run \
    # --publish-dir . \
    # --sample-id $sample \
    # --assembly $assembly \
    # --organism.genus $organism \
    # --reference $ref \
    # --stub

    touch outputs.json
    touch plasmidfinder.json
    mkdir logs
    touch logs/messages.log
    """
}
