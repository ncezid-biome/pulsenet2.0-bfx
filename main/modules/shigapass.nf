process SHIGAPASS {   
    publishDir "$params.publish_dir/${sample}/SHIGAPASS/", pattern: "ShigaPass_summary.csv", mode : "copy"

    input:
        tuple val(sample), path(assembly)
        val(shigapass_db)
    output:
        tuple val(sample), path("ShigaPass_summary.csv"), emit: results
        tuple val(sample), path("shigapass_serotype.txt"), emit: serotype

    script:
        """
        zcat $assembly > ${sample}_assembly.fasta
        echo ${sample}_assembly.fasta > file_path.txt

        ShigaPass.sh \
        -l file_path.txt \
        -o ./ \
        -p $shigapass_db \
        -t 1

        shigapass_result=\$(cat ShigaPass_summary.csv | cut -f 8,9 -d ';' | tail -n 1)
        if [[ \$shigapass_result == SF* ]] && [[ -n \$(echo -n \$shigapass_result | cut -f 2 -d ';') ]]; then
            echo SF\$(echo \$shigapass_result | cut -f 2 -d ';') | tr -d '\n' > shigapass_serotype.txt
        else
            echo \$shigapass_result | cut -f 1 -d ';' | tr -d '\n' > shigapass_serotype.txt
        fi
        """

    stub:
        """
        touch ShigaPass_summary.csv
        """
}
