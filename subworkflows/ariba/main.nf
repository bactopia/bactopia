//
// ariba - Gene identification through local assemblies
//
include { ARIBA_GETREF } from '../../modules/ariba/getref/main'
include { ARIBA_RUN } from '../../modules/ariba/run/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_REPORT } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SUMMARY } from '../../modules/csvtk/concat/main'

workflow ARIBA {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_merged_report = Channel.empty()
    ch_merged_summary = Channel.empty()

    // Build database and run Ariba
    ARIBA_GETREF(params.ariba_db)
    ARIBA_RUN(reads, ARIBA_GETREF.out.db)
    ch_versions = ch_versions.mix(ARIBA_RUN.out.versions)

    ARIBA_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:"${params.ariba_db}-report"], report]}.set{ ch_merge_report }
    CSVTK_CONCAT_REPORT(ch_merge_report, 'tsv', 'tsv')
    ch_merged_report = ch_merged_report.mix(CSVTK_CONCAT_REPORT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_REPORT.out.versions)

    ARIBA_RUN.out.summary.collect{_meta, summary -> summary}.map{ summary -> [[id:"${params.ariba_db}-summary"], summary]}.set{ ch_merge_summary }
    CSVTK_CONCAT_SUMMARY(ch_merge_summary, 'csv', 'csv')
    ch_merged_summary = ch_merged_summary.mix(CSVTK_CONCAT_SUMMARY.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SUMMARY.out.versions)

    emit:
    results = ARIBA_RUN.out.results
    merged_report = ch_merged_report
    merged_summary = ch_merged_summary
    logs = ARIBA_RUN.out.logs.mix(
        CSVTK_CONCAT_REPORT.out.logs,
        CSVTK_CONCAT_SUMMARY.out.logs
    )
    nf_logs = ARIBA_RUN.out.nf_begin.mix(
        ARIBA_RUN.out.nf_err,
        ARIBA_RUN.out.nf_log,
        ARIBA_RUN.out.nf_out,
        ARIBA_RUN.out.nf_run,
        ARIBA_RUN.out.nf_sh,
        ARIBA_RUN.out.nf_trace,
        CSVTK_CONCAT_REPORT.out.nf_begin,
        CSVTK_CONCAT_REPORT.out.nf_err,
        CSVTK_CONCAT_REPORT.out.nf_log,
        CSVTK_CONCAT_REPORT.out.nf_out,
        CSVTK_CONCAT_REPORT.out.nf_run,
        CSVTK_CONCAT_REPORT.out.nf_sh,
        CSVTK_CONCAT_REPORT.out.nf_trace,
        CSVTK_CONCAT_SUMMARY.out.nf_begin,
        CSVTK_CONCAT_SUMMARY.out.nf_err,
        CSVTK_CONCAT_SUMMARY.out.nf_log,
        CSVTK_CONCAT_SUMMARY.out.nf_out,
        CSVTK_CONCAT_SUMMARY.out.nf_run,
        CSVTK_CONCAT_SUMMARY.out.nf_sh,
        CSVTK_CONCAT_SUMMARY.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
