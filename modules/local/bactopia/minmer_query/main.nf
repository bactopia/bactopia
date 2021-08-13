nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "minmer_query"

process MINMER_QUERY {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${sample} - ${dataset_basename}"
    label "max_cpus"
    label "minmer_query"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.txt"

    input:
    tuple val(sample), val(single_end), path(fq), path(sourmash)
    each path(dataset)

    output:
    path "*.txt"
    path "${PROCESS_NAME}/*" optional true

    shell:
    dataset_name = dataset.getName()
    dataset_basename = dataset.getSimpleName()
    mash_w = params.screen_w ? "-w" : ""
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    '''
    LOG_DIR="!{PROCESS_NAME}/!{dataset_basename}"
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
            check-staging.py --fq1 !{fq[0]} --extra !{sourmash} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{sourmash}
        fi
    fi

    if [ "!{dataset_name}" == "refseq-k21-s1000.msh" ]; then
        echo "# Mash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        mash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1

        printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > !{sample}-refseq-k21.txt
        gzip -cd !{fastq} | \
        mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
        sort -gr >> !{sample}-refseq-k21.txt 2> ${LOG_DIR}/mash.err
    elif [ "!{dataset_name}" == "plsdb.msh" ]; then
        echo "# Mash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        mash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1

        printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > !{sample}-plsdb-k21.txt
        gzip -cd !{fastq} | \
        mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset}  - | \
        sort -gr >> !{sample}-plsdb-k21.txt 2> ${LOG_DIR}/mash.err
    elif [ "!{dataset_name}" == "genbank-k21.json.gz" ]; then
        echo "# Sourmash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        sourmash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        sourmash lca classify --query !{sourmash} --db !{dataset} > !{sample}-genbank-k21.txt 2> ${LOG_DIR}/sourmash.err
    elif [ "!{dataset_name}" == "genbank-k31.json.gz" ]; then
        echo "# Sourmash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        sourmash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        sourmash lca classify --query !{sourmash} --db !{dataset} > !{sample}-genbank-k31.txt 2> ${LOG_DIR}/sourmash.err
    else
        echo "# Sourmash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        sourmash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        sourmash lca classify --query !{sourmash} --db !{dataset}  > !{sample}-genbank-k51.txt 2> ${LOG_DIR}/sourmash.err
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
    dataset_name = dataset.getName()
    """
    mkdir ${PROCESS_NAME}
    touch ${sample}.txt
    touch ${PROCESS_NAME}/${sample}
    """
}
