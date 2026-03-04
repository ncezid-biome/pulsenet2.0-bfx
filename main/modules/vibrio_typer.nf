process VIBRIO_SEROTYPER {
   
    input:
        tuple val(sample), path(assembly)
    output:
        tuple val(sample), path("vibrio_serotyper_pcr.json"), path("vibrio_serotyper_results.json"), emit: vibrio_serotyper_jsons
        tuple val(sample), path("outputs.json"), path("master_input.json"), emit: master
        path ("*"), emit: allOut
        tuple val(sample), path("logs/messages.log")

    script:
        """
        ngs-run \
        --sample-id $sample \
        --publish-dir ./ \
        --assembly $assembly\
        --organism.genus VIBRIO.SEROTYPE

        echo "[\$(cat vibrio_serotyper_pcr.json),\$(cat vibrio_serotyper_results.json)]" > master_input.json
        """

    stub:
        """
        touch vibrio_serotyper_pcr.json
        touch vibrio_serotyper_results.json
        touch master_input.json
        touch outputs.json
        mkdir logs
        touch logs/messages.log
        """
}

process VIBRIO_VIRULENCE_FINDER {

    input:
        tuple val(sample), path(assembly)
    output:
        tuple val(sample), path("vibrio_virulence_pcr.json"), path("vibrio_virulence_results.json"), emit: vibrio_virulence_jsons
        tuple val(sample), path("outputs.json"), path("master_input.json"), emit: master
        path ("*"), emit: allOut
        tuple val(sample), path("logs/messages.log")

    script:
        """
        ngs-run \
        --sample-id $sample \
        --publish-dir ./ \
        --assembly $assembly\
        --organism.genus VIBRIO.VIRULENCE
        echo "[\$(cat vibrio_virulence_pcr.json),\$(cat vibrio_virulence_results.json)]" > master_input.json
        """

    stub:
        """
        touch vibrio_virulence_pcr.json
        touch vibrio_virulence_results.json
        touch master_input.json
        touch outputs.json
        mkdir logs
        touch logs/messages.log
        """
}