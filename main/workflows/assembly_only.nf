//load modules
include {ALLELE_CALLING}             from '../modules/AlleleCalling.nf'
include {ALLELE_NAMING}              from '../modules/allele_naming.nf'
include {ANI}                        from '../modules/ani'
include {FETCH_ASSEMBLY_PATH}        from '../modules/fetch_paths'

//load utility functions
include {filterQC}                      from '../utils/utils'
include {getSimilarity}                 from '../utils/utils'
include {GZIP_ASSEMBLY}                 from '../utils/utils'

process INIT_ASSEMBLY_ONLY_MASTER {
    publishDir "$params.publish_dir/${sample}/", pattern: "PipelineProcessOutputs.json", mode : "copy"
    publishDir "$params.publish_dir/${sample}/", pattern: "ADX_analysis_data.json", mode : "copy"
    
    input:
        tuple val(sample), path(assembly)
        path(masterOutScript)

    output:
        tuple val(sample), path(assembly), path("PipelineProcessOutputs.json"), emit: master
        tuple val(sample), path("ADX_analysis_data.json")

    script:
    """
        python3 $masterOutScript --process AssemblyInit
        echo "[]" > ADX_analysis_data.json
    """
}

process ADX_ASSEMBLY_ONLY{
    publishDir "$params.publish_dir/${sample}/", mode : "copy"

    input:
        tuple val(sample), path(alleleNaming)

    output:
        tuple val(sample), path("ADX_analysis_data.json")

    script:
    """
    #!/usr/bin/env python3

    import json
    import gzip

    all_outputs = []

    alleleNaming = json.load(gzip.open("calls_standard.json.gz"))
    allele_out = {}
    allele_out['processType'] = {"name" : "wgmlst"}
    if len(alleleNaming["values"]) > 0:
        allele_out['result'] = {"type": "characterization", "values": alleleNaming["values"] }
    else:
        allele_out['result'] = {"type": "characterization", "values": "Error retrieving alleles" }
    all_outputs.append(allele_out)

    with open('ADX_analysis_data.json', 'w') as output_file:
        json.dump(all_outputs, output_file, indent=4)
    """
}

workflow ASSEMBLY_ONLY {
    if (params.source=="custom") {
        inputs = Channel
                .fromPath(params.reads)
                .splitCsv(header:true, strip:true)
                .map { row -> tuple( row.sample, row.file_path) }
    }
    FETCH_ASSEMBLY_PATH(inputs) 
    GZIP_ASSEMBLY(FETCH_ASSEMBLY_PATH.out.assemblyPath)
    INIT_ASSEMBLY_ONLY_MASTER(GZIP_ASSEMBLY.out, params.masterOutScript)
    assemblies = INIT_ASSEMBLY_ONLY_MASTER.out.master
    
    ANI(params.qckb, assemblies, params.ref_ani, "ANI", params.masterOutScript)
    filtered_ani_orgs = filterQC(ANI.out.organism).map{ tuple(it[0], it[1].trim(), it[2])}
    
    allele_call_input = getSimilarity(assemblies.join(filtered_ani_orgs).map{ tuple(it[0], it[1], it[3], it[4]) })
    ALLELE_CALLING(allele_call_input, params.blast_kb, params.blastdb, params.loci, params.qckb, params.masterOutScript)
    ADX_ASSEMBLY_ONLY(ALLELE_CALLING.out.standard_calls)

}
