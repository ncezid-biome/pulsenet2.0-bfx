//load modules
include {ALLELE_CALLING}             from '../modules/AlleleCalling.nf'
include {ALLELE_FILTERING}           from '../modules/AlleleCalling.nf'
include {ALLELE_NAMING}              from '../modules/allele_naming.nf'
include {MLST}                       from '../modules/mlst'
include {FIND_GENES}                 from '../modules/find_genes'
include {ORG_ANI}                    from '../modules/ani'
include {ISPCR}                      from '../modules/ispcr'
include {SEQSERO}                    from '../modules/seqsero'
include {PLASMIDFINDER}              from '../modules/plasmidfinder'
include {COLLECT_RESULTS_SALM}       from '../modules/collect_results'
include {ADX_ANALYSIS_DATA}          from '../modules/adx_analysis_data'

//load utility functions
include {filterQC}                  from '../utils/utils'
include {removeFromEnd}             from '../utils/utils'
include {compareGenera}             from '../utils/utils'
include {getSimilarity}             from '../utils/utils'

workflow SALM_WORKFLOW {
    
    // [sample, assembly, "SALM", MasterOutput, CRAM, read1, read2, "Salmonella", fastp results, narst_assembly]
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

    //plasmidfinder_input = [sampleID, assembly(narst_assembly), genus(ani)]
    plasmidfinder_input = narst_assembly
        .join(ani_genus)
        .join(qc_end)

    //assembly_ani_genotyping = [sampleID, assembly(corrAssem), genus(ani)]
    assembly_ani_genotyping = assembly
        .join(ani_genus)
        .join(qc_end)
    
    MLST(corrAssem_fullGenus, params.qckb, params.scheme_mapping_kb, params.mlstCC_script, params.cc_mapping)

    corrAssem_fullGenusSpecies = narst_assembly
        .join(midas_genus)
        .join(qc_end)
        .join(master_output)

    FIND_GENES(corrAssem_fullGenusSpecies)

    PLASMIDFINDER(assembly_ani_genotyping, params.plasmid_ref)

    // find_genes_input = [sampleID, assembly[corrAssem], Genus(midas)]
    //find_genes_input = assembly.join(midas_genus)
    //FIND_GENES(find_genes_input)
    
    ISPCR(assembly_ani_genotyping)

    ORG_ANI(params.qckb, assembly_ani_genotyping, params.salm_ani, "ORG_ANI")
    
    
    //seqsero_input = [sampleID, read1, read2, assembly(corrAssem), 
    //                  genus(ani), insilicopcr.json, best_hit.json(org_ani)]
    seqsero_input = reads
        .join(assembly)
        .join(ani_genus)
        .join(ISPCR.out.isPCRjson)
        .join(filterQC(ORG_ANI.out.best)) //QC filter for org ani
        .join(MLST.out.outputs)

    SEQSERO(seqsero_input)


    //Collect most current master to attach genotyping results (currently ani, but will change to allele)
    //latest_master = [sampleID, tempMasterOut(alleleNaming)]
    //latest_master = qcFiltered_assembly_ani.map { tuple( it[0], it[3]) }
    latest_master = ALLELE_FILTERING.out.master.join(qc_end)

    //genotype_master_input = [sampleID, latestMaster(ani), outputs.json(org_ani), best_hit.json(org_ani), 
    //        outputs.json(ispcr), insilicopcr.json(ispcr), outputs.json(seqsero), 
    //        serotype.json(seqsero), outputs.json(findgenes), AMRFinderResultsRaw.json(findgenes),
    //        outputs.json(PF), plasmidfinder.json(PF), outputs.json(MLST)]
    genotype_master_input = latest_master
        .join(ORG_ANI.out.outputs.join(filterQC(ORG_ANI.out.best)))
        .join(ISPCR.out.master)
        .join(SEQSERO.out.master)
        .join(FIND_GENES.out.master)
        .join(PLASMIDFINDER.out.master)
        .join(MLST.out.outputs) // used to print mlst results to master out

    //COLLECT_GENOTYPE_RESULTS(genotype_master_input, params.masterOutScript, params.genotypeFormatter)
    COLLECT_RESULTS_SALM(genotype_master_input, params.masterOutScript, params.genotypeFormatter)

    //adx_input = [masterOut(collect_Results), genotypeResults(collect_Results), fastpReport(prepReads), outputs.json(alleleNaming)]
    adx_input = COLLECT_RESULTS_SALM.out.adx_input
        .join(fastp)
        .join(removeFromEnd(ALLELE_FILTERING.out.standard_calls, 2))
    ADX_ANALYSIS_DATA(adx_input, params.adxScript)
}