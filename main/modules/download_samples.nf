process DOWNLOAD_SAMPLES {

    publishDir "$params.publish_dir/${sample_id}/", pattern: "PipelineProcessOutputs.json", mode : "copy"
    publishDir "$params.publish_dir/${sample_id}/", pattern: "ADX_analysis_data.json", mode : "copy"

    input:
        val(sample_id)
        path(initialMaster)

    output:
        tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz"), path("PipelineProcessOutputs.json"), env(QC), emit: fastq
    script:
        """
        fasterq-dump $sample_id
        gzip *.fastq
        
        QC="true"
        if [[ -e "${sample_id}_1.fastq.gz" ]] && [[ -e "${sample_id}_2.fastq.gz" ]]; then 
            QC="true"
            cp ${initialMaster} PipelineProcessOutputs.json
        else
            QC="false"
            touch ${sample_id}_1.fastq.gz
            touch ${sample_id}_2.fastq.gz
            jq -n '{"sectionName": "DEINTERLEAVE","paired_end": "FALSE", "description": "Sample does not have paired end reads. Paired end reads are required to run the workflow.", "rules":[]}' > issues.json
            jq '.Deinterleave.result="FAIL"' ${initialMaster} > PipelineProcessOutputs_1.json
            jq '.Deinterleave.issues += [input]' PipelineProcessOutputs_1.json issues.json > PipelineProcessOutputs_2.json    
            jq '.metadata.GENUS=""' PipelineProcessOutputs_2.json > PipelineProcessOutputs.json  
        fi
        """
    stub:
        """
        touch raw_deinterleave_1.fastq.gz
        touch raw_deinterleave_2.fastq.gz
        """
}