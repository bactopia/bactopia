//
// ectyper - In-silico prediction of Escherichia coli serotype
//
include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow ECTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ECTYPER_MODULE(fasta)

    // Merge results
    ECTYPER_MODULE.out.tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'ectyper'], tsv]}.set{ ch_merge_ectyper }
    CSVTK_CONCAT(ch_merge_ectyper, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = ECTYPER_MODULE.out.tsv
    txt = ECTYPER_MODULE.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results = ECTYPER_MODULE.out.tsv.mix(
        ECTYPER_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    )
    logs = ECTYPER_MODULE.out.logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = CSVTK_CONCAT.out.nf_begin.mix(
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        ECTYPER_MODULE.out.nf_begin,
        ECTYPER_MODULE.out.nf_err,
        ECTYPER_MODULE.out.nf_log,
        ECTYPER_MODULE.out.nf_out,
        ECTYPER_MODULE.out.nf_run,
        ECTYPER_MODULE.out.nf_sh,
        ECTYPER_MODULE.out.nf_trace
    )
    versions = ECTYPER_MODULE.out.versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
