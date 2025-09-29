//
// tblastx - Search against translated nucleotide BLAST databases using translated nucleotide queries
//
include { BLAST_TBLASTX as TBLASTX_MODULE } from '../../modules/blast/tblastx/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow TBLASTX {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]
    query

    main:
    TBLASTX_MODULE(fasta, query)

    // Merge results
    TBLASTX_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'tblastx'], tsv]}.set{ ch_merge_tblastx }
    CSVTK_CONCAT(ch_merge_tblastx, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = TBLASTX_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = TBLASTX_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
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
    versions = TBLASTX_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
