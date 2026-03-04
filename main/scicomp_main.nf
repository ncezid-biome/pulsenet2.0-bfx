nextflow.enable.dsl=2

//load modules
include {DOWNLOAD_SAMPLES}              from './modules/download_samples'
include {MIDAS as GENUS_IDENTIFY}       from './modules/midas'   
include {PREPARE_READS}                 from './modules/prepare_reads' 
include {GENERATE_ASSEMBLY}             from './modules/generate_assembly'
include {CORRECT_ASSEMBLY}              from './modules/correct_assembly'
include {MIDAS}                         from './modules/midas'
include {ANI}                           from './modules/ani'

//load utility functions
include {filterQC}                      from './utils/utils'
include {removeFromEnd}                 from './utils/utils'
include {compareGenera}                 from './utils/utils'
include {getSimilarity}                 from './utils/utils'

//load subworkflows
include {SALM_WORKFLOW}                 from './workflows/salm'
include {STEC_WORKFLOW}                 from './workflows/stec'
include {LISTERIA_WORKFLOW}             from './workflows/listeria'
include {CAMPY_WORKFLOW}                from './workflows/campy'
include {CRONO_WORKFLOW}                from './workflows/crono'
include {VIBRIO_WORKFLOW}               from './workflows/vibrio'
include {CBOT_WORKFLOW}                 from './workflows/cbot'
include {YERSINIA_WORKFLOW}             from './workflows/yersinia'

//load entrypoint workflows
include {ASSEMBLY_ONLY as assembly}                         from './workflows/assembly_only'
include {STEC_ASSEMBLY_TYPING as stec_assembly_typing}      from './workflows/stec_assembly_typing'

workflow {
    if (params.source=="custom") {
        fastq_path_temp = Channel
                .fromPath(params.reads)
                .splitCsv(header:true, strip:true)
                .map { row -> tuple( row.sample, row.file1, row.file2, params.initialMaster) }
    } 
    else if (params.source == "ncbi"){
        accesssions = Channel.from(file(params.reads).text.tokenize()[1..-1])
        fastq_path_temp = DOWNLOAD_SAMPLES(accesssions, params.initialMaster)
        fastq_path_temp = filterQC(DOWNLOAD_SAMPLES.out.fastq)
    }
    else { 
        println("Invalid input source specified")
    }

    fastq_path = fastq_path_temp.map { tuple (it[0], it[1], it[2], it[3], "DEFAULT") }
    GENUS_IDENTIFY(fastq_path, params.midas_dir, params.qckb, "GENUS_IDENTIFICATION", params.masterOutScript)

    // raw_reads = [sampleID, read1, read2, genus (genus_identify), tempMasterOut (genus_identify)]
    raw_reads = removeFromEnd(fastq_path, 2).join(removeFromEnd(GENUS_IDENTIFY.out.genus, 1))
    PREPARE_READS(raw_reads, params.genome_size, params.qckb, params.masterOutScript)
    
    qcFiltered_reads = filterQC(PREPARE_READS.out.reads).join(removeFromEnd(GENUS_IDENTIFY.out.genus, 2))
    MIDAS(qcFiltered_reads, params.midas_dir, params.qckb, "MIDAS", params.masterOutScript)

    // **midasGenusCheck compares Genus_Identify and Midas genera and filters midas QC**
    // midasGenusCheck = [sampleID, genus(midas), tempMasterOut (midas)]
    midasGenusCheck = compareGenera(removeFromEnd(GENUS_IDENTIFY.out.genus, 2).join(filterQC(MIDAS.out.genus)), "midas")

    // clean_reads = [sampleID, read1, read2, genus(midas), tempMasterOut (midas)]
    cleaned_reads = removeFromEnd(qcFiltered_reads, 2).join(midasGenusCheck)
    
    GENERATE_ASSEMBLY(cleaned_reads, params.qckb, params.masterOutScript)
   
    // reads_assembly = [sampleID, read1, read2, genus (midas), assembly (genAssem), tempMasterOut (genAssem)]
    reads_assembly = removeFromEnd(cleaned_reads, 1).join(filterQC(GENERATE_ASSEMBLY.out.assembly))
    CORRECT_ASSEMBLY(reads_assembly, params.cleaningkb, params.qckb, params.masterOutScript)

    qcFiltered_assembly = filterQC(CORRECT_ASSEMBLY.out.assembly)
    ANI(params.qckb, qcFiltered_assembly, params.ref_ani, "ANI", params.masterOutScript)

    // **midas_ani_genusCheck compares midas and ani and filters ANI QC**
    // midas_ani_genusCheck [sampleID, genus (ani), tempMasterOut (ani)]
    midas_ani_genusCheck = compareGenera(removeFromEnd(midasGenusCheck, 1).join(filterQC(ANI.out.organism)), "ani")

    //qcFiltered_assembly_ani = [sampleID, assembly(corrAssem), genus(ani), tempMasterOut(ani)] (**Input for next process to collect master out**)
    qcFiltered_assembly_ani = removeFromEnd(qcFiltered_assembly, 1).join(midas_ani_genusCheck)
    

    // Separate organisms for subworkflows
    // [sampleID, assembly(corrAssem), genus(ani), genus(midas) tempMasterOut(ani)]
    salm_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "SALM" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)

    stec_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "STEC" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)
        .join(removeFromEnd(raw_reads,2))
        
    listeria_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "LISTERIA" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)
    
    campy_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "CAMPY" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)

    crono_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "CRONO" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)

    vibrio_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "VIBRIO" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)

    cbot_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "CBOT" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)
    
    yersinia_samples = qcFiltered_assembly_ani
        .join(CORRECT_ASSEMBLY.out.alignment)
        .filter{ it[2] == "YERSINIA" }
        .join(removeFromEnd(cleaned_reads, 2))
        .join(removeFromEnd(midasGenusCheck, 1))
        .join(PREPARE_READS.out.fastp)
        .join(CORRECT_ASSEMBLY.out.narst_assembly)

    SALM_WORKFLOW(salm_samples)
    STEC_WORKFLOW(stec_samples)
    LISTERIA_WORKFLOW(listeria_samples)
    CAMPY_WORKFLOW(campy_samples)
    CRONO_WORKFLOW(crono_samples)
    VIBRIO_WORKFLOW(vibrio_samples)
    CBOT_WORKFLOW(cbot_samples)
    YERSINIA_WORKFLOW(yersinia_samples)

}
