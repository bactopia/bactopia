//
// tblastn - Search against translated nucleotide BLAST databases using protein queries
//
include { BLAST_TBLASTN as TBLASTN_MODULE } from '../../modules/blast/tblastn/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow TBLASTN {
    take:
    reads // channel: [ val(meta), [ fasta ] ]
    query

    main:
    ch_versions = Channel.empty()

    // Run TBLASTN
    TBLASTN_MODULE(reads, query)
    ch_versions = ch_versions.mix(TBLASTN_MODULE.out.versions)

    // Merge results
    TBLASTN_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'tblastn'], tsv]}.set{ ch_merge_tblastn }
    CSVTK_CONCAT(ch_merge_tblastn, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    emit:
    tsv = TBLASTN_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = TBLASTN_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = TBLASTN_MODULE.out.nf_begin.mix(
        TBLASTN_MODULE.out.nf_err,
        TBLASTN_MODULE.out.nf_log,
        TBLASTN_MODULE.out.nf_out,
        TBLASTN_MODULE.out.nf_run,
        TBLASTN_MODULE.out.nf_sh,
        TBLASTN_MODULE.out.nf_trace,
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
