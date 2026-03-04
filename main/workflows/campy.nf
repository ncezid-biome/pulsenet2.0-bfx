//load modules
include {ALLELE_CALLING}             from '../modules/AlleleCalling.nf'
include {ALLELE_FILTERING}           from '../modules/AlleleCalling.nf'
include {ALLELE_NAMING}              from '../modules/allele_naming.nf'
include {MLST}                       from '../modules/mlst'
include {FIND_GENES}                 from '../modules/find_genes'
include {PLASMIDFINDER}              from '../modules/plasmidfinder'
include {COLLECT_RESULTS_CAMPY}      from '../modules/collect_results'
include {ADX_ANALYSIS_DATA}          from '../modules/adx_analysis_data'


//load utility functions
include {filterQC}                  from '../utils/utils'
include {removeFromEnd}             from '../utils/utils'
include {compareGenera}             from '../utils/utils'
include {getSimilarity}             from '../utils/utils'

workflow CAMPY_WORKFLOW {
    
    // [sample, assembly, "CAMPY", MasterOutput, CRAM, read1, read2, "Campylobacter", fastp results, narst_assembly]
    take: data

    main:

    // break up data
    reads = data.map{ tuple(it[0], it[5], it[6])}
    assembly = data.map{ tuple( it[0], it[1])}
    cram = data.map{ tuple( it[0], it[4])}
    ani_genus = data.map{ tuple( it[0], it[2])}
    midas_genus = data.map{ tuple( it[0], it[7])}
    master_output = data.map{ tuple( it[0], it[3])}
    fastp = data.map{ tuple( it[0], it[8])}
    narst_assembly = data.map{ tuple( it[0], it[9])}

    //allele_call_input = sampleID, assembly(corrAssem), genus(ani), tempMasterOut(ani), similarity]
    allele_call_input = getSimilarity(assembly.join(ani_genus).join(master_output))
    ALLELE_CALLING(allele_call_input, params.blast_kb, params.blastdb, params.loci, params.qckb, params.masterOutScript)

    //allele_filter_input = [sampleID, assembly(corrAssem), genus(ani), cram_alignment(corrAssem), allele_calls_bam (alleleCall), tempMasterOut(alleleCall)]
    allele_filter_input = assembly
        .join(ani_genus)
        .join(cram)
        .join(filterQC(ALLELE_CALLING.out.allele_calls_bam))
    ALLELE_FILTERING(allele_filter_input, params.filtering_kb, params.qckb, params.masterOutScript)

    //allele_naming_input = [sampleID, allele_calls_xml(alleleCall), standard_calls(alleleFilter), masterOut(alleleFilter)]
    //allele_naming_input = ani_genus.join(ALLELE_CALLING.out.allele_calls_xml).join(filterQC(ALLELE_FILTERING.out.standard_calls))


    //ALLELE_NAMING(allele_naming_input, params.allele_cache_kb, params.allele_cache_ref, params.nomenclature_setting, params.masterOutScript)

    // QC modules end after Allele Filtering and the rest of the modules can now run
    // This is done so non QC modules have to wait until after QC to run
    // in case of a failure

    //qc_end = [sampleID]
    qc_end = removeFromEnd(
        filterQC(ALLELE_FILTERING.out.standard_calls), 
        2
        )

    // corrAssem_fullGenus = [sampleID, assembly[corrAssem], Genus(midas)]
    corrAssem_fullGenus = assembly
        .join(midas_genus)
        .join(qc_end)

    //assembly_ani_genotyping = [sampleID, assembly(corrAssem), genus(ani)]
    assembly_ani_genotyping = narst_assembly
        .join(ani_genus)
        .join(qc_end)
    
    MLST(corrAssem_fullGenus, params.qckb, params.scheme_mapping_kb, params.mlstCC_script, params.cc_mapping)

    corrAssem_fullGenusSpecies = narst_assembly
        .join(midas_genus)
        .join(qc_end)
        .join(master_output)

    FIND_GENES(corrAssem_fullGenusSpecies)

    PLASMIDFINDER(assembly_ani_genotyping, params.plasmid_ref)
    

    //Collect most current master to attach genotyping results
    //latest_master = [sampleID, tempMasterOut(alleleNaming)]
    latest_master = ALLELE_FILTERING.out.master.join(qc_end)

    //genotype_master_input = [sampleID, latestMaster(ani), outputs.json(org_ani), best_hit.json(org_ani), 
    //        outputs.json(ispcr), insilicopcr.json(ispcr), outputs.json(seqsero), 
    //        serotype.json(seqsero), outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes),
    //        outputs.json(PF), plasmidfinder.json(PF), outputs.json(MLST)]
    genotype_master_input = latest_master
        .join(FIND_GENES.out.master)
        .join(PLASMIDFINDER.out.master)
        .join(MLST.out.outputs)

    //COLLECT_GENOTYPE_RESULTS(genotype_master_input, params.masterOutScript, params.genotypeFormatter)
    COLLECT_RESULTS_CAMPY(genotype_master_input, params.masterOutScript, params.genotypeFormatter)

    //adx_input = [masterOut(collect_Results), genotypeResults(collect_Results), fastpReport(prepReads), outputs.json(alleleNaming)]
    adx_input = COLLECT_RESULTS_CAMPY.out.adx_input.join(fastp).join(removeFromEnd(ALLELE_FILTERING.out.standard_calls, 2))
    ADX_ANALYSIS_DATA(adx_input, params.adxScript)
}