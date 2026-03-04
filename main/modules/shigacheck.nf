process SHIGACHECK {
   
    input:
        tuple val(sample), path(isPCR_json), path(PathotypeFinder_output_json), path(VirulenceFinder_json)
        path(shigacheck1_script)
    output:
        tuple val(sample), stdout, emit: result

    script:
        """
        python3 $shigacheck1_script --ispcr $isPCR_json --pf $PathotypeFinder_output_json --vf $VirulenceFinder_json
        """

    stub:
        """
        echo shigella
        """
}
