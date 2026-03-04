process ISPCR {
   
    input:
        tuple val(sample), path(assembly), val(org)
    output:
        tuple val(sample), path("insilicopcr.json"), emit: isPCRjson
        tuple val(sample), path("outputs.json"), path("insilicopcr.json"), emit: master
        path ("*"), emit: allOut
        tuple val(sample), path("logs/messages.log")

    script:
        """
        ngs-run \
        --sample-id $sample \
        --publish-dir ./ \
        --assembly $assembly\
        --organism.genus $org\
        """

    stub:
        """
        # ngs-run \
        # --sample-id $sample \
        # --publish-dir ./ \
        # --assembly $assembly\
        # --organism.genus $org\
        # --stub

        touch insilicopcr.json
        touch outputs.json
        mkdir logs
        touch logs/messages.log
        """
}
