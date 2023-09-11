//
// tblastx - Search against translated nucleotide BLAST databases using translated nucleotide queries
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'tblastx')
options.args = [
    "-outfmt '6 ${params.tblastx_outfmt}'",
    "-qcov_hsp_perc ${params.tblastx_qcov_hsp_perc}",
    "-max_target_seqs ${params.tblastx_max_target_seqs}",
    params.tblastx_opts
].join(' ').replaceAll("\\s{2,}", " ").trim()
QUERY = params.tblastx_query ? file(params.tblastx_query) : []

include { BLAST_TBLASTX as TBLASTX_MODULE } from '../../../modules/nf-core/blast/tblastx/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'tblastx-concat', process_name: params.merge_folder] )

workflow TBLASTX {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run TBLASTX
    TBLASTX_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(TBLASTX_MODULE.out.versions)

    // Merge results
    TBLASTX_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'tblastx'], tsv]}.set{ ch_merge_tblastx }
    CSVTK_CONCAT(ch_merge_tblastx, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = TBLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
