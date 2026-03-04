process SHIGEIFINDER {
     publishDir "$params.publish_dir/${sample}/SHIGEIFINDER/", pattern: "*_result.tsv", mode : "copy"

    input:
        tuple val(sample), path(reads1), path (reads2)
    output:
        tuple val(sample), path("*_result.tsv"), emit: results
        tuple val(sample), path("shigeifinder_serotype.txt"), emit: serotype

    script:
        """
        shigeifinder \
        -t 1 \
        --output ${sample}_result.tsv \
        -r \
        -i $reads1 $reads2

        cat ${sample}_result.tsv | cut -f 5 | tail -n 1 | tr -d '\n' > shigeifinder_serotype.txt
        """

    stub:
        """
        touch ${sample}_result.tsv
        """
}
