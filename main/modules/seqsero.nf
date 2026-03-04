process SEQSERO {

    input:
        tuple val(sample), path(read_one), path(read_two), path(assembly), val(organism), path(isPCR_out), path(ani_out), path(mlst_out, stageAs: "mlst_outputs.json")

    output:
        tuple val(sample), path("outputs.json"), path("serotype.json"), emit: master
        tuple val(sample), path("*"), emit: seqsero2_out
        tuple val(sample), path("logs/messages.log")

    script:
        """
        mkdir reads
        mv $read_one $read_two reads/ # seqsero2s links the reads to ./ so they can't be in ./ to begin with.

        ngs-run \
        --sample-id $sample \
        --publish-dir . \
        --read1 reads/$read_one \
        --read2 reads/$read_two \
        --assembly $assembly \
        --ispcr-out $isPCR_out \
        --species-ani-out $ani_out \
        --organism.genus $organism \
        --mlst-out $mlst_out \
        """

    stub:
        """
        # ngs-run \
        # --sample-id $sample \
        # --publish-dir . \
        # --read1 $read_one \
        # --read2 $read_two \
        # --assembly $assembly \
        # --ispcr-out $isPCR_out \
        # --species-ani-out $ani_out \
        # --organism.genus $organism \
        # --stub

        touch outputs.json
        touch serotype.json
        mkdir logs
        touch logs/messages.log
        """
}

