//load modules
include {SHIGAPASS}                 from '../modules/shigapass.nf'
include {SHIGEIFINDER}              from '../modules/shigeifinder.nf'

//load utility functions
include {removeFromEnd}             from '../utils/utils'

workflow SHIGELLA_WORKFLOW {
    
    // [sample, assembly, read1, read2, True]
    take: data

    main:

    // break up data
    reads = data.map{ tuple(it[0], it[2], it[3])}
    assembly = data.map{ tuple( it[0], it[1])}


    SHIGAPASS(assembly, params.shigapass_db)

    shigeifinder_input = reads
        .join(SHIGAPASS.out.serotype)
        .filter{ it[3].text ==~ /S[BDF]/ || it[3].text =~ /EIEC/ || it[3].text =~ /Unknown/ } // SS doesn't produce serotypes, Unknown is an edge case on SF and possibly others
    

    SHIGEIFINDER {
        removeFromEnd(shigeifinder_input,1)
    }


    emit: 
    shigapass = SHIGAPASS.out.serotype
    shigeifinder = SHIGEIFINDER.out.serotype
}


