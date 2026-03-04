process FETCH_READ_PATHS {

    input:
        tuple val(sample), path(reads), path(tempMaster)

    output:
        tuple val(sample), path("*{_R1_*,_R1,_1}.{fastq,fq}.gz"), path("*{_R2_*,_R2,_2}.{fastq,fq}.gz"), path(tempMaster), emit: fastqPath


        
    script:
    """
    files=(\$(ls $reads))
    num_files="\${#files[@]}"

    if [[ num_files -ne 2 ]]; then
        echo "unexepected number of read files found. Expected 2, but found \$num_files." >&2
        exit 1
    fi

    file1=\${files[0]}
    file2=\${files[1]}

    if [[ $reads == */ ]]; then
        read1="${reads}\$file1"
        read2="${reads}\$file2"
    else
        read1="${reads}/\$file1"
        read2="${reads}/\$file2"
    fi

    ln -s \$read1 .
    ln -s \$read2 .

    """
}

process FETCH_ASSEMBLY_PATH {

    input:
        tuple val(sample), path(assembly_path)

    output:
        tuple val(sample), path("*.{fasta,fa,fna}{,.gz}"), emit: assemblyPath

        
    script:
    """
    file="\$(ls $assembly_path | grep -E '*\\.(fasta|fna|fa)(|\\.gz)')"

    if [[ $assembly_path == */ ]]; then
        assembly="${assembly_path}\$file"
    else
        assembly="${assembly_path}/\$file"
    fi

    ln -s \$assembly .
    """
}

process FETCH_ASSEMBLY_PATH_W_NARST {

    input:
        tuple val(sample), path(file_paths)

    output:
        tuple val(sample), path("*.{fasta,fa,fna}.gz"), emit: assemblyPath
        tuple val(sample), path("NARST_/*.{fasta,fa,fna}.gz"), emit: narstAssemblyPath
        tuple val(sample), path("*{_R1_*,_R1,_1}.{fastq,fq}.gz"), path("*{_R2_*,_R2,_2}.{fastq,fq}.gz"), emit: fastqPath

    script:
    """
    #Grab both assembly files
    file1="\$(ls $file_paths | grep -E '*\\.(fasta|fna|fa)\\.gz' | grep -v 'narst')"
    file2="\$(ls $file_paths | grep -E '*\\.(fasta|fna|fa)\\.gz' | grep 'narst')"

    file3="\$(ls $file_paths | grep -E '*\\.(fastq|fq)\\.gz' | grep -E '(_R1_|_R1\\.|_1\\.)')"
    file4="\$(ls $file_paths | grep -E '*\\.(fastq|fq)\\.gz' | grep -E '(_R2_|_R2\\.|_2\\.)')"

    if [[ $file_paths == */ ]]; then
        assembly="${file_paths}\$file1"
        narst="${file_paths}\$file2"
        read1="${file_paths}\$file3"
        read2="${file_paths}\$file4"
    else
        assembly="${file_paths}/\$file1"
        narst="${file_paths}/\$file2"
        read1="${file_paths}/\$file3"
        read2="${file_paths}/\$file4"
    fi

    ln -s \$read1 .
    ln -s \$read2 .

    cp \$assembly .
    mkdir NARST_
    cp \$narst NARST_/

    """
}
