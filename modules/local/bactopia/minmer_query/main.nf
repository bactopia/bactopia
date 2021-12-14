nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'minmer_query')

process MINMER_QUERY {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${meta.id} - ${dataset_basename}"
    label "base_mem_8gb"
    label "minmer_query"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    input:
    tuple val(meta), path(fq), path(sourmash)
    each path(dataset)

    output:
    path "${meta.id}-${program}-${database}-${kmer}.txt", emit: result
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    dataset_basename = dataset.getSimpleName()
    dataset_info = dataset_basename.split("-") // mash-refseq-k21 --> mash, refseq, k21
    program = dataset_info[0]
    database = dataset_info[1]
    kmer = dataset_info[2]
    mash_w = params.no_winner_take_all ? "" : "-w"
    fastq = meta.single_end ? "${fq[0]}" : "${fq[0]} ${fq[1]}"
    '''
    OUTPUT="!{meta.id}-!{program}-!{database}-!{kmer}.txt"
    OUTPUT_ERR="!{program}-!{database}-!{kmer}.stderr.txt"
    if [ "!{program}" == "mash" ]; then
        echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/\t/g' > ${OUTPUT}
        gzip -cd !{fastq} | \
            mash screen !{mash_w} -i !{params.screen_i} -p !{task.cpus} !{dataset} - | \
            sort -gr >> ${OUTPUT} 2> ${OUTPUT_ERR}
    elif [ "!{program}" == "sourmash" ]; then
        sourmash lca classify --query !{sourmash} --db !{dataset} > ${OUTPUT} 2> ${OUTPUT_ERR}
    fi

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    minmer_query:
        mash: $(echo $(mash --version 2>&1))
        sourmash: $(echo $(sourmash --version 2>&1) | sed 's/sourmash //;')
    END_VERSIONS
    '''
}
