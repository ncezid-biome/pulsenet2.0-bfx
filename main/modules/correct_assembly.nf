process CORRECT_ASSEMBLY {

    publishDir "$params.publish_dir/${sample}/CORRECT_ASSEMBLY/", pattern: "corrected_assembly.fasta", mode : "copy"
    publishDir "$params.publish_dir/${sample}/CORRECT_ASSEMBLY/", pattern: "narst_corrected_assembly.fasta", mode : "copy"
    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"

    input:
        tuple val(sample), path(read1), path(read2), val(org_genus), path(assembly), path(tempMasterOut, stageAs: "tempMasterOut.json")
        path(cleaning_kb)
        path(qc_kb)
        path(masterOutScript)

    output:
        tuple val(sample), path("corrected_assembly.fasta.gz"), path("PipelineProcessOutputs.json"), env(QC), emit: assembly
        tuple val(sample), path("corrected_alignment.cram"), emit: alignment
        tuple val(sample), path("depth_contigs.tsv"), emit: stats
        tuple val(sample), path("IDs.txt"), emit: narst_assembly_ids
        tuple val(sample), path("narst_corrected_assembly.fasta.gz"), emit: narst_assembly
        tuple val(sample), path("outputs.json"), emit: outputs
        tuple val(sample), path("logs/messages.log")
        tuple val(sample), path("corrected_assembly.fasta"), path("narst_corrected_assembly.fasta"), emit: just_for_publish

    script:
        """
        cp $assembly raw_assembly_for_narst.fasta.gz

        ngs-run \
        --sample-id $sample \
        --publish-dir . \
        --n-threads 8 \
        --assembly $assembly \
        --read1 $read1 \
        --read2 $read2 \
        --cleaning-kb.path $cleaning_kb \
        --qc-kb.path $qc_kb \
        --organism.genus $org_genus

        QC="true"
        qcResult="\$(grep -Eo 'PASS|WARN|FAIL|N/A' outputs.json | head -n1)"
        if [[ $params.ignore_qc == "true" && "\$qcResult" != "N/A" ]]; then
            QC="true"
        else
            if [[ "\$qcResult" == "FAIL" || "\$qcResult" == "N/A" ]]; then
                QC="false"
            fi
        fi
        
        awk '\$4 > 10.0 {print \$1}' depth_contigs.tsv > IDs.txt
        seqkit grep -n -f IDs.txt raw_assembly_for_narst.fasta.gz -o narst_corrected_assembly.fasta.gz

        gunzip -k corrected_assembly.fasta.gz
        gunzip -k narst_corrected_assembly.fasta.gz

        python3 $masterOutScript --process CorrectAssembly --tempMaster $tempMasterOut --parseFile outputs.json
        """

    stub:
        """
        # ngs-run \
        # --sample-id $sample \
        # --publish-dir . \
        # --n-threads 8 \
        # --assembly $assembly \
        # --read1 $read1 \
        # --read2 $read2 \
        # --cleaning-kb.path $cleaning_kb \
        # --qc-kb.path $qc_kb \
        # --organism.genus $org_genus
        # --stub

        touch corrected_assembly.fasta.gz
        touch corrected_alignment.cram
        touch depth_contigs.tsv
        touch outputs.json
        touch IDs.txt
        touch narst_corrected_assembly.fasta.gz
            mkdir logs
    touch logs/messages.log

        QC="true"
        touch PipelineProcessOutputs.json
        """
}


