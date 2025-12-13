/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of [Abricate](https://github.com/tseemann/abricate)
 * to screen genome assemblies for antimicrobial resistance and virulence genes, followed by
 * aggregating the results into a single summary report.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, virulence, workflow
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules abricate_run, abricate_summary
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Tab-delimited report of hits per sample
 * @output merged_tsv  Merged report containing hits from all samples
 * @output results     Aggregated channel containing all result files
 * @output logs        Aggregated channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { ABRICATE_RUN     } from '../../modules/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/abricate/summary/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow ABRICATE {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    ABRICATE_RUN(assembly)
    ABRICATE_SUMMARY(gather(ABRICATE_RUN.out.report, 'abricate'))

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = ABRICATE_RUN.out.report
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = ABRICATE_SUMMARY.out.report

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
