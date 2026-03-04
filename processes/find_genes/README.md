# Getting started

This repository includes the following process : 

- FindGenes

To use this process in your pipeline, you must have built and pushed the Docker image to an external registry.
Once the image is available in your registry, you can use it in your nextflow process with the following script.

```groovy
process findGenes {
    container '<your_image_repository>:tag'

    memory { 8.GB * task.attempt }
    cpus { 4 * task.attempt }
    maxRetries 2

    publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/findGenes", mode: params.output_mode_in_use

    input:
        tuple val(id), val(sample_id), path(assembly)
        path(amrfinder_kb)

    output:
        tuple val(id), path("outputs.json"), emit: outputs
        tuple val(id), path("find_genes.tsv"), emit: result

    script:
    """
    ngs-run \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --amrfinder-kb.path $amrfinder_kb
    """

    stub:
    """
    ngs-run \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --amrfinder-kb.path $amrfinder_kb \
    --stub
    """
}  
```

The amrfinder-kb folder should contain the amrfinder.tar.gz knowledge base.

# Running locally
The following examples assume you have the required files locally.  
You can either call using the local source code:
```bash
poetry run ngs-run --sample-id 1 --publish-dir /tmp/ --assembly tests/files/assembly.fasta.gz  --amrfinder-kb.path tests/files/amrfinder_kb --stub
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app ngs-pipeline-process-find-genes ngs-run --sample-id <sample_id> --publish-dir /tmp/ --assembly <assembly.fasta.gz> --amrfinder-kb.path <amrfinder_kb_folder>
```

When running locally you have to use the `--stub` flag as you are likely to be missing some dependencies.  
With the Docker image you can run both in stub and non-stub modes.