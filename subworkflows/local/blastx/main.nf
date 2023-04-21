//
// blastx - Search against protein BLAST databases using translated nucleotide queries
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'blastx')
options.args = [
    "-outfmt '6 ${params.blastx_outfmt}'",
    "-qcov_hsp_perc ${params.blastx_qcov_hsp_perc}",
    "-max_target_seqs ${params.blastx_max_target_seqs}",
    params.blastx_opts
].join(' ').replaceAll("\\s{2,}", " ").trim()
QUERY = params.blastx_query ? file(params.blastx_query) : []

include { BLAST_BLASTX as BLASTX_MODULE } from '../../../modules/nf-core/blast/blastx/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'blastx-concat', process_name: params.merge_folder] )

workflow BLASTX {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run BLASTX
    BLASTX_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(BLASTX_MODULE.out.versions)

    // Merge results
    BLASTX_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'blastx'], tsv]}.set{ ch_merge_blastx }
    CSVTK_CONCAT(ch_merge_blastx, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = BLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
