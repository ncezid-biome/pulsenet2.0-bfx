process DEINTERLEAVE {

  publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"
  publishDir "$params.publish_dir/${sample}/", pattern: "ADX_analysis_data.json", mode : "copy"

  input:
    tuple val(sample), path(read_inter), path(tempMasterOut)

  output:
    tuple val(sample), path("*_1.fastq.gz"), path("*_2.fastq.gz"), path("PipelineProcessOutputs.json"), env(QC), emit: reads
    tuple val(sample), path("ADX_analysis_data.json")
  script:
    """
    echo "[]" > ADX_analysis_data.json
    
    fasterq-dump $read_inter
    gzip *.fastq

        
    QC="true"
    if [[ -e "${read_inter}_1.fastq.gz" ]] && [[ -e "${read_inter}_2.fastq.gz" ]]; then 
        QC="true"
        cp ${tempMasterOut} PipelineProcessOutputs.json
    else
        QC="false"
        touch ${read_inter}_1.fastq.gz
        touch ${read_inter}_2.fastq.gz
        jq -n '{"sectionName": "DEINTERLEAVE","paired_end": "FALSE", "description": "Sample does not have paired end reads. Paired end reads are required to run the workflow.", "rules":[]}' > issues.json
        jq '.Deinterleave.result="FAIL"' ${tempMasterOut} > PipelineProcessOutputs_1.json
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