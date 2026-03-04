# Getting started

This repository includes one algorithm : 

- GenerateAssembly

To use this algorithm in your pipeline's processes, you must have built & pushed Docker image to an external registry.
Once the image available in your registry, you can use it in your nextflow process with the following script.

**GenerateAssembly**

```groovy
process generateAssembly {
  container '<your_image_repository>:tag'

    memory { 8.GB * task.attempt }
    cpus { 4 * task.attempt }
    errorStrategy { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    maxRetries 2

    publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/generateAssembly", mode: params.output_mode_in_use

    input:
        tuple val(id), val(sample_id), val(genus), val(species), path(read1), path(read2)
        path(qc_kb)

    output:
        tuple val(id), path("outputs.json"), emit: outputs
        tuple val(id), path("assembly.fasta.gz"), emit: assembly

    script:
    """
    ngs-run \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --read1 $read1 \
    --read2 $read2 \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus}
    """

    stub:
    """
    ngs-run \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --read1 $read1 \
    --read2 $read2 \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus} \
    --stub
    """
}
``` 

# Running locally
The following examples assume you have the required files locally.  
You can either call using the local source code:
```bash
poetry run ngs-run --sample-id 1 --publish-dir /tmp/ --read1 tests/files/read1.fastq --read2 tests/files/read2.fastq --qc-kb.path tests/files/qc_kb/ --organism.genus SALMONELLA --stub
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app ngs-pipeline-process-generate-assembly ngs-run --sample-id <sample_id> --publish-dir /tmp/ --read1 <read1.fastq.gz> --read2 <read2.fastq.gz> --qc-kb.path <qc_kb_folder> --organism.genus <genus> --organism.species <species>
```

When running locally you have to use the `--stub` flag as you are likely to be missing some dependencies.  
With the Docker image you can run both in stub and non-stub modes.