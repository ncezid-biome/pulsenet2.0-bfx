//load modules
include {FETCH_ASSEMBLY_PATH_W_NARST}   from '../modules/fetch_paths'
include {ANI}                           from '../modules/ani'
include {MLST}                          from '../modules/mlst'
include {FIND_GENES}                    from '../modules/find_genes'
include {ISPCR}                         from '../modules/ispcr'
include {PLASMIDFINDER}                 from '../modules/plasmidfinder'
include {VIRULENCE_FINDER}              from '../modules/virulencefinder'
include {KMA}                           from '../modules/kma'
include {PATHOTYPE_FINDER}              from '../modules/pathotypefinder'
include {COLLECT_RESULTS_STEC}          from '../modules/collect_results'
include {STXTYPER}                      from '../modules/stxtyper.nf'
include {STXCONDENSER}                  from '../modules/stxcondenser.nf'
include {ADX_ANALYSIS_DATA}             from '../modules/adx_analysis_data'

//load utility functions
include {filterQC}                      from '../utils/utils'
include {removeFromEnd}                 from '../utils/utils'

process INIT_STEC_ASSEMBLY_MASTER {
    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"
    
    input:
        tuple val(sample), path(assembly)
        path(masterOutScript)

    output:
        tuple val(sample), path(assembly), path("PipelineProcessOutputs.json"), emit: master

    script:
    """
        python3 $masterOutScript --process stecassemblytyping
    """
}

process ADX_STEC_ASSEMBLY_TYPING_ONLY {
    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(genotypingResults)

    output:
        tuple val(sample), path("ADX_analysis_data.json")
    
    script:
    """
    #!/usr/bin/env python3

    import json
    import gzip
    import re

    all_outputs = []

    genotype_input = json.load(open("GenotypingResult.json"))
    genotype_content = open("GenotypingResult.json").read()

    ## Plasmid Results
    plasmid_out = {}
    plasmid_out['processType'] = {"name" : "plasmids"}
    oldVal = genotype_input['plasmids.json']['results']
    newVal = { key: 1 if value else 2 for key, value in oldVal.items() }
    plasmid_out['result'] = {"type": "characterization", "values": newVal}
    all_outputs.append(plasmid_out)

    ## Resistance (AMRFinder/Find Genes) Results
    resist_out = {}
    resist_out['processType'] = {"name" : "resistance"}
    oldVal = genotype_input['resistance.json']['results']
    newVal = { key: 1 if value else 2 for key, value in oldVal.items() }
    resist_out['result'] = {"type": "characterization", "values": newVal}
    all_outputs.append(resist_out)

    ## GenotypingResults Full File
    jsonString = re.sub(r"\\s+", " ", genotype_content).replace("{ ", "{")
    geno_out = {}
    geno_out['processType'] = {"name" : "genotypingResult"}
    geno_out['result'] = {"type": "genotypingResult", "rawJsonAnalysis": jsonString}
    all_outputs.append(geno_out)

    with open('ADX_analysis_data.json', 'w') as output_file:
        json.dump(all_outputs, output_file, indent=4)
    """
}

workflow STEC_ASSEMBLY_TYPING {
    if (params.source=="custom") {
        inputs = Channel
                .fromPath(params.reads)
                .splitCsv(header:true, strip:true)
                .map { row -> tuple( row.sample, row.file_path) }
    }
    FETCH_ASSEMBLY_PATH_W_NARST(inputs)
    //FETCH_ASSEMBLY_PATH_W_NARST.out.assemblyPath.view()
    //FETCH_ASSEMBLY_PATH_W_NARST.out.narstAssemblyPath.view()
    //FETCH_ASSEMBLY_PATH_W_NARST.out.fastqPath.view()

    //master = INIT_STEC_ASSEMBLY_MASTER(FETCH_ASSEMBLY_PATH_W_NARST.out.assemblyPath, params.masterOutScript)
    //    .map { tuple( it[0], it[2]) }

    //master.view()
    
    INIT_STEC_ASSEMBLY_MASTER(FETCH_ASSEMBLY_PATH_W_NARST.out.assemblyPath, params.masterOutScript)
    
    ANI(params.qckb, INIT_STEC_ASSEMBLY_MASTER.out.master, params.ref_ani, "ANI", params.masterOutScript)
    
    
    //[val(sample), stdout, path("PipelineProcessOutputs.json"), env(QC)]
    filtered_ani_orgs = filterQC(ANI.out.organism).map{ tuple(it[0], it[1].trim(), it[2])}
    
    //collect latest master and qc end
    latest_master = filtered_ani_orgs.map { tuple(it[0], it[2]) }
    qc_end = removeFromEnd(filtered_ani_orgs, 2)

    //set data
    assembly = FETCH_ASSEMBLY_PATH_W_NARST.out.assemblyPath.join(qc_end)
    narst_assembly = FETCH_ASSEMBLY_PATH_W_NARST.out.narstAssemblyPath.join(qc_end)
    cleaned_reads = FETCH_ASSEMBLY_PATH_W_NARST.out.fastqPath.join(qc_end)
    ani_genus = filtered_ani_orgs.map { tuple(it[0], it[1]) }
    midas_genus = qc_end.map { tuple(it[0], "ESCHERICHIA") }
    
    
    // corrAssem_fullGenus = [sampleID, assembly[corrAssem], Genus(midas)]
    corrAssem_fullGenus = assembly
        .join(midas_genus)

    findgenes_amr_input = narst_assembly
        .join(midas_genus)
        .join(latest_master)

    plasmidfinder_input = narst_assembly
        .join(ani_genus)

    //assembly_ani_genotyping = [sampleID, assembly(corrAssem), genus(ani)]
    assembly_ani_genotyping = assembly
        .join(ani_genus)

    MLST(corrAssem_fullGenus, params.qckb, params.scheme_mapping_kb, params.mlstCC_script, params.cc_mapping)

    FIND_GENES(findgenes_amr_input)

    PLASMIDFINDER(plasmidfinder_input, params.plasmid_ref)
    VIRULENCE_FINDER(assembly_ani_genotyping, params.virulence_ref)
    
    ISPCR(assembly_ani_genotyping)
    
    kma_refs = tuple([params.stec_otypes, params.stec_htypes, params.stec_lookup_table])
    KMA(cleaned_reads, kma_refs)

    PATHOTYPE_FINDER(assembly_ani_genotyping, params.pathotypefinder_db)

    STXTYPER(cleaned_reads, params.holotoxins)

    STXCONDENSER(STXTYPER.out.stxtypejson.join(ISPCR.out.isPCRjson))

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

    COLLECT_RESULTS_STEC(genotype_master_input, params.masterOutScript, params.genotypeFormatter)

    //adx_input = [sample, genotypeResults(collect_Results)]
    adx_input = COLLECT_RESULTS_STEC.out.adx_input.map { tuple(it[0], it[2]) }
    //ADX_STEC_ASSEMBLY_TYPING_ONLY(adx_input, "/scicomp/home-pure/lsz0/PN2.0/NewStructure/main/pn2.0_test/testGenotyping.py")
    ADX_STEC_ASSEMBLY_TYPING_ONLY(adx_input)
}
