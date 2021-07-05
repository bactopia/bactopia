nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "download_references"

process DOWNLOAD_REFERENCES {
    /*
    Download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${sample} - ${params.max_references} reference(s)"
    label "download_references"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: 'mash-dist.txt'

    input:
    tuple val(sample), val(single_end), path(fq), path(sample_sketch)
    path(refseq_sketch)

    output:
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("genbank/*.gbk"), emit:CALL_VARIANTS_AUTO, optional: true
    path("mash-dist.txt")
    file "${PROCESS_NAME}/*" optional true

    shell:
    no_cache = params.no_cache ? '-N' : ''
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    total = params.max_references
    '''
    LOG_DIR="!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    # Verify AWS files were staged
    if [[ ! -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "true" ]; then
            check-staging.py --fq1 !{fq[0]} --extra !{sample_sketch} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{sample_sketch}
        fi
    fi

    # Get Mash distance
    echo "# Mash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    mash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    mash dist -t !{sample_sketch} !{refseq_sketch} | grep -v "query" | sort -k 2,2 > distances.txt

    # Pick genomes to download
    printf "accession\\tdistance\\tlatest_accession\\tupdated\\n" > mash-dist.txt
    select-references.py distances.txt !{total} !{tie_break} >> mash-dist.txt

    # Pick only latest accessions
    grep -v distance mash-dist.txt | cut -f3 > download-list.txt

    # Download genomes
    echo "# ncbi-genome-download Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    ncbi-genome-download --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    ncbi-genome-download bacteria -l complete -o ./ -F genbank -p !{task.cpus} -A download-list.txt -r !{params.max_retry} !{no_cache} > ${LOG_DIR}/ncbi-genome-download.out 2> ${LOG_DIR}/ncbi-genome-download.err

    # Move and uncompress genomes
    mkdir genbank_temp
    find refseq -name "*.gbff.gz" | xargs -I {} mv {} genbank_temp/
    rename 's/(GC[AF]_\\d+).*/$1/' genbank_temp/*
    mkdir genbank
    ls genbank_temp/ | xargs -I {} sh -c 'gzip -cd genbank_temp/{} > genbank/!{sample}-{}.gbk'
    rm -rf genbank_temp

    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate GenBank files
        rm -rf refseq/
    fi

    # pass the FASTQs along
    mkdir -p fastqs
    if [[ -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        fi
    else
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            cp !{fq[0]} fastqs/!{sample}_R1.fastq.gz
            cp !{fq[1]} fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            cp  !{fq[0]} fastqs/!{sample}.fastq.gz
        fi
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{PROCESS_NAME}.err
        cp .command.out ${LOG_DIR}/!{PROCESS_NAME}.out
        cp .command.sh ${LOG_DIR}/!{PROCESS_NAME}.sh || :
        cp .command.trace ${LOG_DIR}/!{PROCESS_NAME}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir fastqs
    mkdir genbank
    mkdir ${PROCESS_NAME}
    touch fastqs/${sample}.fastq.gz
    touch genbank/*.gbk
    touch ${PROCESS_NAME}/${sample}
    touch mash-dist.txt
    """
}
