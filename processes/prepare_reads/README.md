# Getting started

This repository includes two algorithms : 

- Download Reads
- Cleanup Reads

To use these algorithms in your pipeline's processes, you must have built & pushed Docker image to an external registry.
Once the image available in your registry, you can use it in your nextflow process with the following scripts.

Note that the parameter `aws` for downloadReads is a boolean and as such needs to be passed as a flag rather than a key/value pair.
For that we use a ternary operator based on its value, it is present only when it evaluates to true.  
Likewise the parameter `species` for cleanupReads is only passed if it is defined (in stub mode it can be omitted altogether).

**Download Reads**

```groovy
process downloadReads {
  container '<your_image_repository>:tag'

  publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/downloadReads", mode: params.output_mode_in_use

  input:
    tuple val(id), val(sample_id), val(accession_id)
    val(aws)
  
  output:
    tuple val(id), path("outputs.json"), emit: outputs
    tuple val(id), path("read_1.fastq.gz"), path("read_2.fastq.gz"), emit: reads

  script:
    """
    ngs-run DownloadReads \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --accession-id $accession_id \
    ${aws ? '--aws': ''}
     """

  stub:
    """
    ngs-run DownloadReads \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --accession-id $accession_id \
    ${aws ? '--aws': ''} \
    --stub
    """
}
```

**Cleanup Reads**
```groovy
process cleanupReads {
  container '<your_image_repository>:tag'


  publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/cleanupReads", mode: params.output_mode_in_use


  input:
    tuple val(id), val(sample_id), path(read1), path(read2), val(genus), val(species)
    path(genome_kb)
    path(qc_kb)
  
  output:
    tuple val(id), path("outputs.json"), emit: outputs
    tuple val(id), path("fastp_report.json"), emit: report
    tuple val(id), path("cleaned_read_1.fastq.gz"), path("cleaned_read_2.fastq.gz"), emit: reads

  script:
    """
    ngs-run CleanupReads \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --read1 $read1 \
    --read2 $read2 \
    --genome-kb.path $genome_kb \
    --qc-kb.path $qc_kb
    """

  stub:
    """
    ngs-run CleanupReads \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --read1 $read1 \
    --read2 $read2 \
    --genome-kb.path $genome_kb \
    --qc-kb.path $qc_kb \
    --stub
    """
}
```

# Running locally
The following examples assume you have the required files locally.  
You can either call using the local source code:
```bash
poetry run ngs-run DownloadReads --stub  --sample-id 1 --publish-dir /tmp/ --aws --accession-id 1
poetry run ngs-run CleanupReads --sample-id 1 --publish-dir /tmp/ --organism.genus fake --read1 tests/files/read1.fastq.gz --read2 tests/files/read2.fastq.gz --genome-kb.path tests/files/genome_kb/ --qc-kb.path tests/files/qc_kb/ --stub
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app ngs-pipeline-process-prepare-reads ngs-run DownloadReads --sample-id <sample_id> --publish-dir /tmp/ --accession-id <accession_id> --aws
docker run -v $(pwd):/app ngs-pipeline-process-prepare-reads ngs-run CleanupReads --sample-id <sample_id> --publish-dir /tmp/ --organism.genus <genus> --read1 <read1.fastq.gz> --read2 <read2.fastq.gz> --genome-kb.path <genome_kb_folder> --qc-kb.path <qc_kb_folder>
```

When running locally you have to use the `--stub` flag as you are likely to be missing some dependencies.  
With the Docker image you can run both in stub and non-stub modes.