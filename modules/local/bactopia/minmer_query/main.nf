nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "minmer_query"

process MINMER_QUERY {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${sample} - ${dataset_basename}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(sample), val(single_end), path(fq), path(sourmash)
    each path(dataset)

    output:
    path "${sample}-${program}-${database}-${kmer}.txt", emit: result
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    dataset_name = dataset.getName()
    dataset_basename = dataset.getSimpleName()
    dataset_info = dataset_basename.split("-") // mash-refseq-k21 --> mash, refseq, k21
    program = dataset_info[0]
    database = dataset_info[1]
    kmer = dataset_info[2]
    mash_w = params.no_winner_take_all ? "" : "-w"
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    '''
    OUTPUT="!{sample}-!{program}-!{database}-!{kmer}.txt"
    OUTPUT_ERR="!{program}-!{database}-!{kmer}.stderr.txt"
    if [ "!{program}" == "mash" ]; then
        printf "identity\tshared-hashes\tmedian-multiplicity\tp-value\tquery-ID\tquery-comment\n" > ${OUTPUT}
        gzip -cd !{fastq} | \
            mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset} - | \
            sort -gr >> ${OUTPUT} 2> ${OUTPUT_ERR}

        # Capture version
        mash --version > mash.version.txt 2>&1
    elif [ "!{program}" == "sourmash" ]; then
        sourmash lca classify --query !{sourmash} --db !{dataset} > ${OUTPUT} 2> ${OUTPUT_ERR}

        # Capture version
        sourmash --version > sourmash.version.txt 2>&1
    fi
    '''

    stub:
    dataset_name = dataset.getName()
    """
    touch ${sample}-mash-refseq-k21.txt
    touch ${sample}-mash-k21.stderr.txt
    """
}
