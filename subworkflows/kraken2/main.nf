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
 * @modules kraken2 as kraken2_module
 *
 * @input reads
 * Channel containing reads data
 *
 * @input database
 * Channel containing database data
 *
 * @output classified     Classified
 * @output kraken2_report Kraken2 Report
 * @output unclassified   Unclassified
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution logs from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { KRAKEN2 as KRAKEN2_MODULE } from '../../modules/kraken2/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow KRAKEN2 {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    KRAKEN2_MODULE(reads, database)

    emit:
    // Individual outputs
    classified: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.classified
    kraken2_report: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.kraken2_report
    unclassified: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.unclassified

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KRAKEN2_MODULE.out.classified,
        KRAKEN2_MODULE.out.kraken2_report,
        KRAKEN2_MODULE.out.unclassified
    ])
    logs: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = KRAKEN2_MODULE.out.versions
}
