//
// clermontyping - in silico phylotyping of Escherichia genus
//
include { CLERMONTYPING as CLERMONTYPING_MODULE } from '../../modules/clermontyping/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow CLERMONTYPING {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    CLERMONTYPING_MODULE(fasta)

    // Merge results
    CLERMONTYPING_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'clermontyping'], tsv]}.set{ ch_merge_clermontyping }
    CSVTK_CONCAT(ch_merge_clermontyping, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = CLERMONTYPING_MODULE.out.tsv
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = CLERMONTYPING_MODULE.out.tsv.mix(
        CSVTK_CONCAT.out.csv
    )
    logs = CLERMONTYPING_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CLERMONTYPING_MODULE.out.nf_begin.mix(
        CLERMONTYPING_MODULE.out.nf_err,
        CLERMONTYPING_MODULE.out.nf_log,
        CLERMONTYPING_MODULE.out.nf_out,
        CLERMONTYPING_MODULE.out.nf_run,
        CLERMONTYPING_MODULE.out.nf_sh,
        CLERMONTYPING_MODULE.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = CLERMONTYPING_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
