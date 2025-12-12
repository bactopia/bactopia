/**
 * Rapidly identify genes by creating local assemblies from paired-end reads.
 *
 * This subworkflow uses [ARIBA](https://github.com/sanger-pathogens/ariba)
 * (Antimicrobial Resistance Identification By Assembly) to rapidly identify genes
 * in a database by creating local assemblies. It first downloads and prepares an ARIBA database,
 * then analyzes paired-end reads to identify genes, and finally aggregates results across all samples.
 *
 * @status stable
 * @keywords bacteria, reads, antimicrobial resistance, virulence, local assembly
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent, resource-download
 * @citation ariba
 *
 * @modules ariba_getref, ariba_run, csvtk_concat
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end reads in FASTQ format
 *
 * @input db
 * Database name for ARIBA analysis (e.g., ncbi, card, vfdb, resfinder, argannot)
 *
 * @output report         Per-sample tab-delimited reports of gene findings
 * @output summary        Per-sample CSV summaries of gene detection results
 * @output merged_report  Consolidated report containing gene findings from all samples
 * @output merged_summary Consolidated summary containing detection results from all samples
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
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
