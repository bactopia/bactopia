//
// ariba - Gene identification through local assemblies
//
include { ARIBA_GETREF } from '../../modules/ariba/getref/main'
include { ARIBA_RUN    } from '../../modules/ariba/run/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_REPORT } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SUMMARY } from '../../modules/csvtk/concat/main'

workflow ARIBA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    db

    main:
    // Build database and run Ariba
    ARIBA_GETREF(db)
    ARIBA_RUN(reads, ARIBA_GETREF.out.db)

    // Merge results
    ARIBA_RUN.out.report.collect{_meta, report -> report}.map{ report -> [[id:"${db}-report", args:'-C "$" --lazy-quotes'], report]}.set{ ch_merge_report }
    CSVTK_CONCAT_REPORT(ch_merge_report, 'tsv', 'tsv')

    ARIBA_RUN.out.summary.collect{_meta, summary -> summary}.map{ summary -> [[id:"${db}-summary", args:'--lazy-quotes'], summary]}.set{ ch_merge_summary }
    CSVTK_CONCAT_SUMMARY(ch_merge_summary, 'csv', 'csv')

    emit:
    // Individual outputs
    report = ARIBA_RUN.out.report
    summary = ARIBA_RUN.out.summary
    merged_report = CSVTK_CONCAT_REPORT.out.csv
    merged_summary = CSVTK_CONCAT_SUMMARY.out.csv

    // Generic aggregate outputs
    results = ARIBA_RUN.out.report.mix(
        ARIBA_RUN.out.summary,
        ARIBA_RUN.out.supplemental,
        CSVTK_CONCAT_REPORT.out.csv,
        CSVTK_CONCAT_SUMMARY.out.csv
    )
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
    versions = ARIBA_RUN.out.versions.mix(
        CSVTK_CONCAT_REPORT.out.versions,
        CSVTK_CONCAT_SUMMARY.out.versions
    )
}
