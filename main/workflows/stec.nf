//load modules
include {ALLELE_CALLING}            from '../modules/AlleleCalling.nf'
include {ALLELE_FILTERING}          from '../modules/AlleleCalling.nf'
include {ALLELE_NAMING}             from '../modules/allele_naming.nf'
include {MLST}                      from '../modules/mlst'
include {FIND_GENES}                from '../modules/find_genes'
include {ISPCR}                     from '../modules/ispcr'
include {PLASMIDFINDER}             from '../modules/plasmidfinder'
include {VIRULENCE_FINDER}          from '../modules/virulencefinder'
include {KMA}                       from '../modules/kma'
include {PATHOTYPE_FINDER}          from '../modules/pathotypefinder'
include {COLLECT_RESULTS_STEC}      from '../modules/collect_results'
include {SHIGACHECK}                from '../modules/shigacheck'
include {STXTYPER}                  from '../modules/stxtyper.nf'
include {STXCONDENSER}              from '../modules/stxcondenser.nf'
include {ADX_ANALYSIS_DATA}         from '../modules/adx_analysis_data'

// load subworkflows
include {SHIGELLA_WORKFLOW}         from './shigella'

//load utility functions
include {filterQC}                  from '../utils/utils'
include {removeFromEnd}             from '../utils/utils'
include {compareGenera}             from '../utils/utils'
include {getSimilarity}             from '../utils/utils'

workflow STEC_WORKFLOW {
    
    // [sample, assembly, "STEC", MasterOutput, CRAM, read1, read2, "Escherichia", narst_assembly, raw_read1, raw_read2 ]
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
    raw_reads = data.map{ tuple(it[0], it[10], it[11])}
    

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
    //allele_naming_input = ALLELE_CALLING.out.allele_calls_xml.join(filterQC(ALLELE_FILTERING.out.standard_calls))
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

    findgenes_amr_input = narst_assembly
        .join(midas_genus)
        .join(qc_end)
        .join(master_output)

    //plasmidfinder_input = [sampleID, assembly(narst_assembly), genus(ani)]
    plasmidfinder_input = narst_assembly
        .join(ani_genus)
        .join(qc_end)

    //assembly_ani_genotyping = [sampleID, assembly(corrAssem), genus(ani)]
    assembly_ani_genotyping = assembly
        .join(ani_genus)
        .join(qc_end)

    MLST(corrAssem_fullGenus, params.qckb, params.scheme_mapping_kb, params.mlstCC_script, params.cc_mapping)

    FIND_GENES(findgenes_amr_input)

    PLASMIDFINDER(plasmidfinder_input, params.plasmid_ref)
    VIRULENCE_FINDER(assembly_ani_genotyping, params.virulence_ref)

    // find_genes_input = [sampleID, assembly[corrAssem], Genus(midas)]
    //find_genes_input = assembly.join(midas_genus)
    //FIND_GENES(find_genes_input)
    
    ISPCR(assembly_ani_genotyping)
    
    kma_refs = tuple([params.stec_otypes, params.stec_htypes, params.stec_lookup_table])
    KMA(reads.join(qc_end), kma_refs)

    PATHOTYPE_FINDER(assembly_ani_genotyping, params.pathotypefinder_db)

    STXTYPER(reads.join(qc_end), params.holotoxins)

    STXCONDENSER(STXTYPER.out.stxtypejson.join(ISPCR.out.isPCRjson))

    
    // shigella check
    shigacheck1_input = ISPCR.out.isPCRjson // sample, ispcrjson
        .join(PATHOTYPE_FINDER.out.master) // sample, outputs_json, pathotypefinder_json
        .join(VIRULENCE_FINDER.out.master) // sample, outputs_json, virulencefinder_json
        .map { tuple(it[0],it[1],it[2],it[5]) } // sample, isPCRjson, outputs_json(pf), VirulenceFinderjson
    SHIGACHECK(shigacheck1_input, params.shigacheck1_script)
   

    shigella_samples = assembly
        .join(raw_reads)
        .join(SHIGACHECK.out.result)
        .filter{ it[4] == "shigella" }

    SHIGELLA_WORKFLOW(shigella_samples)

    //Collect most current master to attach genotyping results
    //latest_master = [sampleID, tempMasterOut(alleleNaming)]
    latest_master = ALLELE_FILTERING.out.master.join(qc_end) //joining with qc_end makes sure to only proceed if QC was passed

    //genotype_master_input = [sampleID, latestMaster(ani), outputs.json(org_ani), best_hit.json(org_ani), 
    //        outputs.json(ispcr), insilicopcr.json(ispcr),
    //        outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes),
    //        outputs.json(PF), plasmidfinder.json(PF), outputs.json(MLST),
    //        outputs.json(PathotypeFinder)]
    genotype_master_input = latest_master
        .join(ISPCR.out.master)
        .join(FIND_GENES.out.master)
        .join(PLASMIDFINDER.out.master)
        .join(MLST.out.outputs)
        .join(PATHOTYPE_FINDER.out.master)
        .join(KMA.out.master)
        .join(STXTYPER.out.master)
        .join(STXCONDENSER.out.master)
        .join(VIRULENCE_FINDER.out.master)
        .join(SHIGELLA_WORKFLOW.out.shigapass, remainder: true)
        .join(SHIGELLA_WORKFLOW.out.shigeifinder, remainder: true)
        .map{ sample, tempMasterOut, ispcr_outputs, insilicopcr, amr_outputs, amr_resultsRaw, pf_outputs, plasmidfinder, mlst_outputs, pathotype_outputs, pathotype_genotypes, kma_outputs, kma_serotype, stx_outputs, stx_results, stx_condenser_outputs, stx_condenser_master, virulence_outputs, virulence_results, shigapass_results, shigeifinder_results -> [sample, tempMasterOut, ispcr_outputs, insilicopcr, amr_outputs, amr_resultsRaw, pf_outputs, plasmidfinder, mlst_outputs, pathotype_outputs, pathotype_genotypes, kma_outputs, kma_serotype, stx_outputs, stx_results, stx_condenser_outputs, stx_condenser_master, virulence_outputs, virulence_results, shigapass_results ?: [], shigeifinder_results ?: [] ] }

    // SHIGELLA workflow must come before COLLECT_RESULTS so that ecoli group field can be written with shigella results
    //COLLECT_GENOTYPE_RESULTS(genotype_master_input, params.masterOutScript, params.genotypeFormatter)
    COLLECT_RESULTS_STEC(genotype_master_input, params.masterOutScript, params.genotypeFormatter)

    //adx_input = [masterOut(collect_Results), genotypeResults(collect_Results), fastpReport(prepReads), outputs.json(alleleNaming)]
    adx_input = COLLECT_RESULTS_STEC.out.adx_input
        .join(fastp)
        .join(removeFromEnd(ALLELE_FILTERING.out.standard_calls, 2))
    ADX_ANALYSIS_DATA(adx_input, params.adxScript)

}


