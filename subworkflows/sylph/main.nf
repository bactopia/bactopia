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
 * @modules sylph_profile
 *
 * @input reads
 * Channel containing reads data
 *
 * @input database
 * Channel containing database data
 *
 * @output tsv      Tsv
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SYLPH_PROFILE } from '../../modules/sylph/profile/main'
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow SYLPH {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    SYLPH_PROFILE(reads, database)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.tsv
    logs: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.versions
}
