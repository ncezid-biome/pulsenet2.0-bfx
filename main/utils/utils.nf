def filterQC(inQC) {
    return inQC
    .filter{it[-1] == "true"}
    .map { it[0..it.size()-2] }
}

def removeFromEnd(inChan, numRemove) {
    return inChan
    .map { it[0..it.size()-(numRemove+1)] }
}

//inChan = [sampleID, genus1, genus2, tempMasterOut]
def compareGenera(inChan, compType) {
    GENUS_ABBREVS = [
        "Salmonella":"SALM",
        "Escherichia":"STEC",
        "Listeria":"LISTERIA",
        "Vibrio":"VIBRIO",
        "Campylobacter":"CAMPY",
        "Cronobacter":"CRONO",
        "Yersinia": "YERSINIA",
        "Clostridium": "CBOT"
    ]
    if (compType == "midas") {
        return inChan
        .filter{ it[1] == it[2] }
        .map { it[0,2,3] }
    }
    else if (compType == "ani") {
        return inChan
        .filter{ GENUS_ABBREVS[it[1]] == it[2].trim() }
        .map { tuple (it[0], it[2].trim(), it[3]) }
    }
    else {
        println "Unknown comparison type: '$compType'."
    }
}

//(sampleID, assembly, organism, tempMaster)
def getSimilarity(inChan) {
    blast_similarity = [
        "SALM":75,
        "STEC":85,
        "LISTERIA":85,
        "VIBRIO":85,
        "CAMPY":70,
        "CRONO":80,
        "CBOT":85,
        "YERSINIA":85
    ]
    return inChan
    .map { tuple (it[0], it[1], it[2], it[3], blast_similarity[it[2]]) }
}

process GZIP_ASSEMBLY {
    stageInMode = 'copy'
    input:
        tuple val(sample), path (file)
    output:
        tuple val(sample), path("assembly.fasta.gz")

    script:
    """
    outfile=assembly.fasta.gz # Add gzip extension if missing
    if [[ $file == *.gz ]];
    then
        # File already gzipped
        gunzip $file -c | tr -d "\\r" > temp.fasta
        gzip temp.fasta -c > \$outfile 
    else
        # Gzipping file
        tr -d "\\r" < $file > temp.fasta
        gzip temp.fasta -c > \$outfile

    fi
    """
}