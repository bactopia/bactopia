/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules abricate_run, abricate_summary
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output tsv        Tsv
 * @output merged_tsv Merged Tsv
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { ABRICATE_RUN     } from '../../modules/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow ABRICATE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    ABRICATE_RUN(fasta)
    ABRICATE_SUMMARY(gather(ABRICATE_RUN.out.report, 'abricate'))

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = ABRICATE_RUN.out.report
    merged_tsv: Channel<Tuple<Map, Path>> = ABRICATE_SUMMARY.out.report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRICATE_RUN.out.report,
        ABRICATE_SUMMARY.out.report
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRICATE_RUN.out.logs,
        ABRICATE_SUMMARY.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRICATE_RUN.out.nf_logs,
        ABRICATE_SUMMARY.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRICATE_RUN.out.versions,
        ABRICATE_SUMMARY.out.versions
    ])
}
