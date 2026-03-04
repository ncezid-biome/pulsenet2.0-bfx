process KMA {
   
    input:
        tuple val(sample), path(read1), path(read2)
        tuple path(o_types), path(h_types), path(lookup)
    output:
        tuple val(sample), path("ec_kma_serotyper_result_out.json"), emit: KMAjson
        tuple val(sample), path("outputs.json"), path("ec_kma_serotyper_result_out.json"), emit: master
        path ("*"), emit: allOut
        tuple val(sample), path("logs/messages.log")

    script:
        """
        ngs-run \
        --sample-id $sample \
        --publish-dir ./ \
        --read1 $read1 \
        --read2 $read2 \
        --otypes $o_types \
        --htypes $h_types \
        --lookuptable $lookup
        """

    stub:
        """
        touch ec_kma_serotyper_result_out.json
        touch outputs.json
        mkdir logs
        touch logs/messages.log
        """
}
