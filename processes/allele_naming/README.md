# Getting started

This repository includes the following process : 

- AlleleNaming

## Use in nextflow

To use these processes in your pipeline, you must have built & pushed Docker image to an external registry.
Once the image available in your registry, you can use it in your nextflow process with the following scripts.

**Name Alleles**

```groovy
process alleleNaming {

    publishDir "${params.output_dir}/${workflow.runName}_${workflow.sessionId}/${(params.use_sample_id) ? sample_id : id}/alleleNaming", mode: params.output_mode_in_use

    input:
        tuple val(hashed_id), val(id), path(cache_kb), path(nomenclature_settings), path(allele_calls_xml), path(standard_calls)

    output:
        tuple val(id), path("outputs.json"), emit: outputs

    script:
    """
    ngs-run \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --allele-calls-xml $allele_calls_xml \
    --allele-calls-profile $standard_calls \
    --allele-cache-kb.path $allele_cache_kb \
    --nomenclature-settings $nomenclature_settings
    """

    stub:
    """
    ngs-run \
    --sample-id $sample_id \
    --publish-dir $task.publishDir.path \
    --allele-calls-xml $allele_calls_xml \
    --allele-calls-profile $standard_calls \
    --allele-cache-kb.path $allele_cache_kb \
    --nomenclature-settings $nomenclature_settings \
    --stub
    """
}
``` 


* The `allele-calls-xml` is the XML output file (gzipped) emitted by the `allele-calling` process.

* The `allele-calls-profile` is the standard calls json file (gzipped) emitted by either the `allele-calling` or the `allele-filtering` process.

* The `allele-cache-kb` is a knowledge base folder, containing a cache folder with at least two files:
  - `accepted_alleles_link`: with content a rotation link (`0` or `1`)
  - `accepted_alleles_{rotation_link}.fasta.gz`: rotation fasta.gz file containing the accepted alleles.

  This cache folder can be copied from the BioNumerics CE shared folder.

* The `nomenclature-settings` is a json file containing the settings to access the BioNumercs CE Nomenclature service:
  ```json
  {
        "url": "<url of the BN CE>",
        "project": "<BN CE project>",
        "password": "<BN CE password>",
        "serial": "<BN serial number that can access the BN CE project>",
        "organism_id": "<organism id used in BN CE>",
        "lab_id": "<lab id that submits new alleles (optional)>"
  }
  ```
  The `organism_id` is typically the prefix given to the wgMLST locus IDs.  
  The `lab_id` is optional, it is added to a comment field in hte nomenclature database when adding a new allele. The default value is `lab_id`.


# Running locally
The following examples assume you have the required files locally.  
You can either call using the local source code, but you have to use the `--stub` flag as you are likely to be missing some dependencies:  
```bash
poetry run ngs-run --sample-id <sample_id> --publish-dir /tmp/ --allele-calls-xml <allele_calls.xml.gz> --allele-calls-profile <calls_standard.json.gz> --allele-cache-kb.path <allele_cache_kb_folder> --nomenclature-settings <settings.json> --stub
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app <docker-image> ngs-run --sample-id <sample_id> --publish-dir /tmp/ --allele-calls-xml <allele_calls.xml.gz> --allele-calls-profile <calls_standard.json.gz> --allele-cache-kb.path <allele_cache_kb_folder> --nomenclature-settings <settings.json>
```
With the Docker image you can run both in stub and non-stub modes.  
The `-v ($pwd):/app` argument means that your local version of the code is going to overwrite the one in the image, so you can quickly test changes.