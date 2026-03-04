process STXTYPER {
   
    input:
        tuple val(sample), path(read1), path(read2)
        path(holotoxins)
    output:
        tuple val(sample), path("ecstxtyper_result_out.json"), emit: stxtypejson
        tuple val(sample), path("outputs.json"), path("ecstxtyper_result_out.json"), emit: master
        path ("*"), emit: allOut
        tuple val(sample), path("logs/messages.log")

    script:
        """
        ngs-run \
        --sample-id $sample \
        --publish-dir ./ \
        --read1 $read1 \
        --read2 $read2 \
        --holotoxins $holotoxins
        """

    stub:
        """
        touch ecstxtyper_result_out.json
        touch outputs.json
        mkdir logs
        touch logs/messages.log
        """
}
