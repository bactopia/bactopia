//
// ectyper - In-silico prediction of Escherichia coli serotype
//
include { ECTYPER as ECTYPER_MODULE } from '../../modules/ectyper/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow ECTYPER {
    take:
    fasta // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ECTYPER_MODULE(fasta)
    ch_versions = ch_versions.mix(ECTYPER_MODULE.out.versions.first()
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions))
    ch_logs = ch_logs.mix(ECTYPER_MODULE.out.logs)

    // Merge results
    ECTYPER_MODULE.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'ectyper'], tsv]}.set{ ch_merge_ectyper }
    CSVTK_CONCAT(ch_merge_ectyper, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)

    emit:
    tsv = ECTYPER_MODULE.out.tsv
    txt = ECTYPER_MODULE.out.txt
    merged_tsv = CSVTK_CONCAT.out.csv
    logs = ch_logs
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
    versions = ch_versions
}
