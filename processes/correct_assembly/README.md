# Getting started

This repository includes one algorithm : 

- CorrectAssembly

To use this algorithm in your pipeline's processes, you must have built & pushed Docker image to an external registry.
Once the image available in your registry, you can use it in your nextflow process with the following script.

**CorrectAssembly**

```groovy
process correctAssembly {
  container '<your_image_repository>:tag'

  memory 8.GB
  cpus 4 

  publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/correctAssembly", mode: params.output_mode_in_use

  input:  
      tuple val(id), val(sample_id), val(organism), path(read1), path(read2), path(assembly)
      path(cleaning_kb)
      path(qc_kb)

  output:
      tuple val(id), path("outputs.json"), emit: outputs
      tuple val(id), path("corrected_assembly.fasta.gz"), emit: assembly
      tuple val(id), path("corrected_alignment.cram"), emit: alignment
      tuple val(id), path("depth_contigs.tsv"), emit: stats

  script:
  """
  ngs-run \
  --sample-id $sample_id \
  --publish-dir $task.publishDir.path \
  --read1 $read1 \
  --read2 $read2 \
  --assembly $assembly \
  --cleaning-kb.path $cleaning_kb \
  --qc-kb.path $qc_kb \
  --organism.genus $organism_genus \
  ${species ? '--organism.species ' + $species : ''}
  """

  stub:
  """
  ngs-run \
  --sample-id $sample_id \
  --publish-dir $task.publishDir.path \
  --read1 $read1 \
  --read2 $read2 \
  --assembly $assembly \
  --cleaning-kb.path $cleaning_kb \
  --qc-kb.path $qc_kb \
  --organism.genus $organism_genus \
  ${species ? '--organism.species ' + $species : ''}
  --stub
  """
}
``` 

# Running locally
The following examples assume you have the required files locally.  
You can either call using the local source code:
```bash
poetry run ngs-run --sample-id 1 --publish-dir /tmp/ --read1 tests/files/read1.fastq.gz --read2 tests/files/read2.fastq.gz --assembly tests/files/assembly.fasta.gz --cleaning-kb.path tests/files/cleaning_kb/ --qc-kb.path tests/files/qc_kb/ --organism.genus fake --stub
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app ngs-pipeline-process-correct-assembly ngs-run --sample-id <sample_id> --publish-dir /tmp/ --read1 <read1.fastq.gz> --read2 <read2.fastq.gz> --assembly <assembly.fasta.gz> --cleaning-kb.path <cleaning_kb_folder> --qc-kb.path <qc_kb_folder> --organism.genus <genus>
```

When running locally you have to use the `--stub` flag as you are likely to be missing some dependencies.  
With the Docker image you can run both in stub and non-stub modes.