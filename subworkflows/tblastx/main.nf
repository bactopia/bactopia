//
// tblastx - Search against translated nucleotide BLAST databases using translated nucleotide queries
//
include { BLAST_TBLASTX as TBLASTX_MODULE } from '../../modules/blast/tblastx/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow TBLASTX {
    take:
    reads // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    QUERY = params.tblastx_query ? file(params.tblastx_query) : []

    // Run TBLASTX
    TBLASTX_MODULE(reads, QUERY)
    ch_versions = ch_versions.mix(TBLASTX_MODULE.out.versions)

    // Merge results
    TBLASTX_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'tblastx'], tsv]}.set{ ch_merge_tblastx }
    CSVTK_CONCAT(ch_merge_tblastx, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = TBLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = TBLASTX_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = TBLASTX_MODULE.out.nf_begin.mix(
        TBLASTX_MODULE.out.nf_err,
        TBLASTX_MODULE.out.nf_log,
        TBLASTX_MODULE.out.nf_out,
        TBLASTX_MODULE.out.nf_run,
        TBLASTX_MODULE.out.nf_sh,
        TBLASTX_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
