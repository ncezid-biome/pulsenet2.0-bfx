//load modules
include {MLST}                          from '../modules/mlst'
include {ALLELE_CALLING}                from '../modules/AlleleCalling.nf'
include {ALLELE_FILTERING}              from '../modules/AlleleCalling.nf'
include {COLLECT_RESULTS_YERSINIA}      from '../modules/collect_results'
include {ADX_ANALYSIS_DATA_QUALITY}     from '../modules/adx_analysis_data'

//load utility functions
include {filterQC}                  from '../utils/utils'
include {removeFromEnd}             from '../utils/utils'
include {compareGenera}             from '../utils/utils'
include {getSimilarity}             from '../utils/utils'




workflow YERSINIA_WORKFLOW {
    
    // [sample, assembly, "Yersinia", MasterOutput, CRAM, read1, read2, "Yersinia", fastp results, narst_assembly]
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

    //qc_end = [sampleID]
    qc_end = removeFromEnd(
        filterQC(ALLELE_FILTERING.out.standard_calls), 
        2
        )

    // corrAssem_fullGenus = [sampleID, assembly[corrAssem], Genus(midas)]
    corrAssem_fullGenus = assembly
        .join(midas_genus)
        .join(qc_end)
    
    // MLST(corrAssem_fullGenus, params.qckb, params.scheme_mapping_kb, params.mlstCC_script, params.cc_mapping)
    
    //Collect most current master to attach genotyping results
    latest_master = ALLELE_FILTERING.out.master.join(qc_end)

    //genotype_master_input = [sampleID, latestMaster(ani), outputs.json(MLST)]
    // genotype_master_input = latest_master
    //     .join(MLST.out.outputs)

    // COLLECT_RESULTS_YERSINIA(genotype_master_input, params.masterOutScript)
    
    // adx_input = COLLECT_RESULTS_YERSINIA.out.adx_input
    adx_input = latest_master
        .join(fastp)
        .join(removeFromEnd(ALLELE_FILTERING.out.standard_calls, 2))

    ADX_ANALYSIS_DATA_QUALITY(adx_input, params.adxScript)
}
