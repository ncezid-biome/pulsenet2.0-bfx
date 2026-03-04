process COLLECT_RESULTS_SALM {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        //[sampleID, latestMaster(ani), outputs.json(org_ani), best_hit.json(org_ani), outputs.json(ispcr), insilicopcr.json(ispcr), outputs.json(seqsero), serotype.json(seqsero), outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes), outputs.json(PF), plasmidfinder.json(PF)]
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(org_ani_outputs, stageAs: "org_ani_outputs.json"),
            path(org_ani_best_hit, stageAs: "org_ani_best_hit.json"),
            path(ispcr_outputs, stageAs: "ispcr_outputs.json"),
            path(insilicopcr),
            path(seqsero_outputs, stageAs: "seqsero_outputs.json"),
            path(serotype),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json")
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process ORG_ANI --tempMaster $tempMasterOut --parseFile $org_ani_outputs --parseAddition $org_ani_best_hit
    python3 $masterOutScript --process AMRFinder --tempMaster PipelineProcessOutputs.json --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process isPCR --tempMaster PipelineProcessOutputs.json --parseFile $ispcr_outputs --parseAddition $insilicopcr
    python3 $masterOutScript --process SeqSero --tempMaster PipelineProcessOutputs.json --parseFile $seqsero_outputs --parseAddition $serotype
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    python3 $masterOutScript --process MLST --tempMaster PipelineProcessOutputs.json --parseFile $mlst_outputs

    python3 $genotypeFormatter \
    --genus salmonella \
    --serotype $serotype \
    --sero-out seqsero_outputs.json \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --ispcr $insilicopcr \
    --ispcr-out ispcr_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json
    

    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}

process COLLECT_RESULTS_STEC {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(ispcr_outputs, stageAs: "ispcr_outputs.json"),
            path(insilicopcr),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json"),
            path(pathotype_outputs, stageAs: "pathotype_outputs.json"),
            path(pathotype_genotypes),
            path(kma_outputs, stageAs: "kma_outputs.json"),
            path(kma_serotype),
            path(stx_outputs, stageAs: "stx_outputs.json"),
            path(stx_results),
            path(stx_condenser_outputs, stageAs: "stx_condenser_outputs.json"),
            path(stx_condenser_master, stageAs: "stx_condenser_master_input.json"),
            path(virulence_outputs, stageAs: "virulence_outputs.json"),
            path(virulence_results),
            path(shigapass_results),
            path(shigeifinder_results)
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    def shigapass_command = shigapass_results ? "python3 $masterOutScript --process shigapass --tempMaster PipelineProcessOutputs.json --parseFile $shigapass_results": ""
    def shigeifinder_command = shigeifinder_results ? "python3 $masterOutScript --process shigeifinder --tempMaster PipelineProcessOutputs.json --parseFile $shigeifinder_results": ""
    """
    python3 $masterOutScript --process AMRFinder --tempMaster $tempMasterOut --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process isPCR --tempMaster PipelineProcessOutputs.json --parseFile $ispcr_outputs --parseAddition $insilicopcr
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    python3 $masterOutScript --process MLST --tempMaster PipelineProcessOutputs.json --parseFile $mlst_outputs
    python3 $masterOutScript --process STXtyper --tempMaster PipelineProcessOutputs.json --parseFile $stx_outputs --parseAddition $stx_results
    python3 $masterOutScript --process STXcondenser --tempMaster PipelineProcessOutputs.json --parseFile $stx_condenser_outputs --parseAddition $stx_condenser_master

    python3 $masterOutScript --process PathotypeFinder --tempMaster PipelineProcessOutputs.json --parseFile $pathotype_outputs
    python3 $masterOutScript --process KMA --tempMaster PipelineProcessOutputs.json --parseFile $kma_outputs --parseAddition $kma_serotype
    python3 $masterOutScript --process virulencefinder --tempMaster PipelineProcessOutputs.json --parseFile $virulence_outputs --parseAddition $virulence_results

    $shigapass_command
    $shigeifinder_command

    #StecGroup needs to occur after KMA and PathotypeFinder
    python3 $masterOutScript --process StecGroup --tempMaster PipelineProcessOutputs.json

    python3 $genotypeFormatter \
    --genus ecoli \
    --serotype $kma_serotype \
    --sero-out kma_outputs.json \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --ispcr $insilicopcr \
    --ispcr-out ispcr_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json \
    --pathotype $pathotype_outputs $pathotype_genotypes \
    --pathotype-out pathotype_outputs.json \
    --stxtyper $stx_results \
    --stx-out stx_outputs.json \
    --stxcondenser $stx_condenser_master \
    --virulence $virulence_results \
    --virulence-out virulence_outputs.json 

    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}

process COLLECT_RESULTS_LISTERIA {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json")
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process AMRFinder --tempMaster $tempMasterOut --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    python3 $masterOutScript --process MLST --tempMaster PipelineProcessOutputs.json --parseFile $mlst_outputs


    python3 $genotypeFormatter \
    --genus listeria \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json
    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}

process COLLECT_RESULTS_CAMPY {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json")
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process AMRFinder --tempMaster $tempMasterOut --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    python3 $masterOutScript --process MLST --tempMaster PipelineProcessOutputs.json --parseFile $mlst_outputs


    python3 $genotypeFormatter \
    --genus campy \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json
    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}

process COLLECT_RESULTS_VIBRIO {
    publishDir "$params.publish_dir/${sample}/", mode : "copy"
    input:
      // [sampleID, latestMaster(allele_naming), outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes),
        //        outputs.json(PF), plasmidfinder.json(PF), outputs.json(MLST), outputs.json(serotype),
        //        master_input.json(serotype), outputs.json(virulence), master_input.json(virulence)]
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json"),
            path(serotype_outputs, stageAs: "serotype_outputs.json"),
            path(serotype_results, stageAs: "serotype_merged.json"),
            path(virulence_outputs, stageAs: "virulence_outputs.json"),
            path(virulence_results, stageAs: "virulence_merged.json")
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process AMRFinder --tempMaster $tempMasterOut --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    python3 $masterOutScript --process MLST --tempMaster PipelineProcessOutputs.json --parseFile $mlst_outputs
    python3 $masterOutScript --process VibrioSerotype --tempMaster PipelineProcessOutputs.json --parseFile serotype_outputs.json --parseAddition serotype_merged.json
    python3 $masterOutScript --process VibrioVirulence --tempMaster PipelineProcessOutputs.json --parseFile virulence_outputs.json --parseAddition virulence_merged.json


    python3 $genotypeFormatter \
    --genus vibrio \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json \
    --serotype $serotype_results \
    --sero-out serotype_outputs.json \
    --virulence $virulence_results \
    --virulence-out virulence_outputs.json
    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}


process COLLECT_RESULTS_CBOT {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
    //[sampleID, latestMaster(ani), outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes), outputs.json(PF), plasmidfinder.json(PF), outputs.json(mlst), outputs.json(cbot toxin typing), toxintyping results]
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json"),
            path(cbot_toxin_outputs, stageAs: "cbot_toxin_outputs.json"),
            path(cbot_toxin_results)
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process cbot_toxintyper --tempMaster $tempMasterOut --parseFile $cbot_toxin_outputs --parseAddition $cbot_toxin_results
    python3 $masterOutScript --process MLST --tempMaster PipelineProcessOutputs.json --parseFile $mlst_outputs
    python3 $masterOutScript --process AMRFinder --tempMaster PipelineProcessOutputs.json --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    
    python3 $genotypeFormatter \
    --genus cbot \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json \
    --ctoxin $cbot_toxin_results \
    --cbot-toxin-out cbot_toxin_outputs.json 

    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}

process COLLECT_RESULTS_CRONO {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
    //[sampleID, latestMaster(ani), outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes), outputs.json(PF), plasmidfinder.json(PF), outputs.json(mlst)]
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(amr_outputs, stageAs: "amr_outputs.json"),
            path(amr_resultsRaw),
            path(pf_outputs, stageAs: "pf_outputs.json"),
            path(plasmidfinder),
            path(mlst_outputs, stageAs: "mlst_outputs.json")
        path(masterOutScript)
        path(genotypeFormatter)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), path("GenotypingResult.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process MLST --tempMaster $tempMasterOut --parseFile $mlst_outputs
    python3 $masterOutScript --process AMRFinder --tempMaster PipelineProcessOutputs.json --parseFile $amr_outputs --parseAddition $amr_resultsRaw
    python3 $masterOutScript --process PlasmidFinder --tempMaster PipelineProcessOutputs.json --parseFile $pf_outputs --parseAddition $plasmidfinder
    
    python3 $genotypeFormatter \
    --genus crono \
    --plasmid $plasmidfinder \
    --plasmid-out pf_outputs.json \
    --amr $amr_resultsRaw \
    --amr-out amr_outputs.json 
    """

    stub:
    """
    touch PipelineProcessOutputs.json
    touch GenotypingResult.json
    """
}

process COLLECT_RESULTS_YERSINIA {

    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
    //[sampleID, latestMaster(ani), outputs.json(mlst)]
        tuple val(sample), path(tempMasterOut, stageAs: "tempMasterOut.json"),
            path(mlst_outputs, stageAs: "mlst_outputs.json")
        path(masterOutScript)

    output:
        tuple val(sample), path("PipelineProcessOutputs.json"), emit: adx_input

    script:
    """
    python3 $masterOutScript --process MLST --tempMaster $tempMasterOut --parseFile $mlst_outputs
    
    """

    stub:
    """
    touch PipelineProcessOutputs.json
    """
}
