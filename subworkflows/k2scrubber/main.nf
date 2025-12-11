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
 * @modules wget, kraken2
 *
 * @input reads
 * Channel containing reads data
 *
 * @output human                Human
 * @output kraken2_report       Kraken2 Report
 * @output scrub_report         Scrub Report
 * @output scrub_special_report Scrub Special Report
 * @output scrubbed             Scrubbed
 * @output scrubbed_extra       Scrubbed Extra
 * @output results              Aggregated results channel containing all output files
 * @output logs                 Aggregated logs channel containing all execution logs
 * @output nf_logs              Aggregated Nextflow execution logs from all processes
 * @output versions             Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { WGET         } from '../../modules/wget/main'
include { KRAKEN2      } from '../../modules/kraken2/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow K2SCRUBBER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    WGET([
        "name": "k2scrubber",
        "save_as": "k2_HPRC_20230810.tar.gz",
        "url": "https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz?download=1"
    ])
    KRAKEN2(reads, WGET.out.download)

    emit:
    // Individual outputs
    human: Channel<Tuple<Map, Path>> = KRAKEN2.out.classified
    kraken2_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.kraken2_report
    scrub_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.scrub_report
    scrub_special_report: Channel<Tuple<Map, Path>> = KRAKEN2.out.scrub_special_report
    scrubbed: Channel<Tuple<Map, Path>> = KRAKEN2.out.unclassified
    scrubbed_extra: Channel<Tuple<Map, Path>> = KRAKEN2.out.unclassified_extra

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        KRAKEN2.out.classified,
        KRAKEN2.out.kraken2_report,
        KRAKEN2.out.scrub_report,
        KRAKEN2.out.unclassified
    ])
    logs: Channel<Tuple<Map, Path>> = KRAKEN2.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = KRAKEN2.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = KRAKEN2.out.versions
}
