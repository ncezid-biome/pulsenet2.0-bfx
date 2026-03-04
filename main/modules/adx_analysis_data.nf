process ADX_ANALYSIS_DATA {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(masterOut), path(genotypeResults), path(fastp), path(alleleFiltering)
        path(adx_script)

    output:
        tuple val(sample), path("ADX_analysis_data.json")

    script:
    """
    python3 $adx_script \
    --masterOutFile $masterOut \
    --genotypeResults $genotypeResults \
    --fastpFile $fastp \
    --wgmlst $alleleFiltering
    """
}

//For organisms that do not run Allele calling or any genotype modules
process ADX_ANALYSIS_DATA_QUALITY {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(masterOut), path(fastp), path(alleleFiltering)
        path(adx_script)

    output:
        tuple val(sample), path("ADX_analysis_data.json")

    script:
    """
    python3 $adx_script \
    --masterOutFile $masterOut \
    --fastpFile $fastp \
    --wgmlst $alleleFiltering
    """
}
