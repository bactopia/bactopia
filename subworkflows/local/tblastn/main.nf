//
// tblastn - Search against translated nucleotide BLAST databases using protein queries
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'tblastn')
options.args = [
    "-outfmt '6 ${params.tblastn_outfmt}'",
    "-qcov_hsp_perc ${params.tblastn_qcov_hsp_perc}",
    "-max_target_seqs ${params.tblastn_max_target_seqs}",
    params.tblastn_opts
].join(' ').replaceAll("\\s{2,}", " ").trim()
QUERY = params.tblastn_query ? file(params.tblastn_query) : []

include { BLAST_TBLASTN as TBLASTN_MODULE } from '../../../modules/nf-core/blast/tblastn/main' addParams( options: options )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'tblastn-concat', process_name: params.merge_folder] )

workflow TBLASTN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    // Run TBLASTN
    TBLASTN_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(TBLASTN_MODULE.out.versions)

    // Merge results
    TBLASTN_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'tblastn'], tsv]}.set{ ch_merge_tblastn }
    CSVTK_CONCAT(ch_merge_tblastn, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = TBLASTN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    versions = ch_versions // channel: [ versions.yml ]
}
