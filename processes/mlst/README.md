# Getting started

This repository includes the following algorithm : 

- MLST

To use this algorithm in your pipeline's processes, you must have built and pushed the Docker image to an external registry.
Once the image is available in your registry, you can use it in your nextflow process with the following script.

```groovy
  process mlst {

    container '<your_image_repository>:tag'

    memory 8.GB
    cpus 2

    publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/mlst", mode: params.output_mode_in_use
    input:
        tuple val(hashed_id), val(id), val(genus), val(species), path(assembly)
        path(qc_kb), stageAs: "qc_kb"
        path(scheme_mapping_kb), stageAs: "scheme_mapping_kb"

    output:
        tuple val(id), path("outputs.json"), emit: outputs
        tuple val(id), path("novel.fasta"), emit: novel_alleles

    script:
    """
    ngs-run \
    --sample-id $id \
    --assembly $assembly \
    --qc-kb.path $qc_kb \
    --scheme-mapping-kb.path $scheme_mapping_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --publish-dir $task.publishDir.path
    """

    stub:
    """
    ngs-run \
    --sample-id $id \
    --assembly $assembly \
    --qc-kb.path $qc_kb \
    --scheme-mapping-kb.path $scheme_mapping_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --publish-dir $task.publishDir.path \
    --stub
    """
}  

# Running locally
The following examples assume you have the required file locally.  
You can either call using the local source code:
```bash
poetry run ngs-run --sample-id test --assembly tests/files/assembly.fasta.gz --stub --publish-dir /tmp/ --qc-kb.path local/knowledge_bases/qc/v1.0.2/ --scheme-mapping-kb.path local/knowledge_bases/scheme-mapping/v1.0.0/  --organism.genus SALMONELLA
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app ngs-pipeline-process-mlst ngs-run --sample-id test --assembly tests/files/assembly.fasta.gz --publish-dir /tmp/  --qc-kb.path local/knowledge_bases/qc/v1.0.2/ --scheme-mapping-kb.path local/knowledge_bases/scheme-mapping/v1.0.0/  --organism.genus SALMONELLA
```

When running locally you have to use the `--stub` flag as you are likely to be missing some dependencies.  
With the Docker image you can run both in stub and non-stub modes.