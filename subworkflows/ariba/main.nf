//
// ariba - Gene identification through local assemblies
//
nextflow.preview.types = true

include { ARIBA_GETREF                         } from '../../modules/ariba/getref/main'
include { ARIBA_RUN                            } from '../../modules/ariba/run/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_REPORT  } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SUMMARY } from '../../modules/csvtk/concat/main'
include { flattenPaths                         } from 'plugin/nf-bactopia'
include { gather                               } from 'plugin/nf-bactopia'

workflow ARIBA {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    db: String

    main:
    ARIBA_GETREF(db)
    ARIBA_RUN(reads, ARIBA_GETREF.out.db)
    CSVTK_CONCAT_REPORT(gather(ARIBA_RUN.out.report, "${db}-report", '-C "$" --lazy-quotes'), 'tsv', 'tsv')
    CSVTK_CONCAT_SUMMARY(gather(ARIBA_RUN.out.summary, "${db}-summary", '--lazy-quotes'), 'csv', 'csv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Path>> = ARIBA_RUN.out.report
    summary: Channel<Tuple<Map, Path>> = ARIBA_RUN.out.summary
    merged_report: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_REPORT.out.csv
    merged_summary: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_SUMMARY.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ARIBA_RUN.out.report,
        ARIBA_RUN.out.summary,
        ARIBA_RUN.out.supplemental,
        CSVTK_CONCAT_REPORT.out.csv,
        CSVTK_CONCAT_SUMMARY.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ARIBA_RUN.out.logs,
        CSVTK_CONCAT_REPORT.out.logs,
        CSVTK_CONCAT_SUMMARY.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ARIBA_RUN.out.nf_logs,
        CSVTK_CONCAT_REPORT.out.nf_logs,
        CSVTK_CONCAT_SUMMARY.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ARIBA_RUN.out.versions,
        CSVTK_CONCAT_REPORT.out.versions,
        CSVTK_CONCAT_SUMMARY.out.versions
    ])
}
