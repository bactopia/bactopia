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
 * @modules srahumanscrubber_initdb, srahumanscrubber_scrub
 *
 * @input reads
 * Channel containing reads data
 *
 * @output scrubbed             Scrubbed
 * @output scrubbed_extra       Scrubbed Extra
 * @output scrub_report         Scrub Report
 * @output scrub_special_report Scrub Special Report
 * @output results              Aggregated results channel containing all output files
 * @output logs                 Aggregated logs channel containing all execution logs
 * @output nf_logs              Aggregated Nextflow execution logs from all processes
 * @output versions             Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB  } from '../../modules/srahumanscrubber/scrub/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SRAHUMANSCRUBBER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)

    emit:
    // Individual outputs
    scrubbed: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    scrubbed_extra: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrubbed_extra
    scrub_report: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    scrub_special_report: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrub_special_report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SRAHUMANSCRUBBER_SCRUB.out.scrubbed,
        SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    ])
    logs: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.versions
}
