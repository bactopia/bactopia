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
 * @modules abritamr_run, csvtk_concat
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output summary_tsv        Summary Tsv
 * @output merged_summary_tsv Merged Summary Tsv
 * @output matches_tsv        Matches Tsv
 * @output partials_tsv       Partials Tsv
 * @output virulence_tsv      Virulence Tsv
 * @output amrfinder_tsv      Amrfinder Tsv
 * @output results            Aggregated results channel containing all output files
 * @output logs               Aggregated logs channel containing all execution logs
 * @output nf_logs            Aggregated Nextflow execution logs from all processes
 * @output versions           Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { ABRITAMR_RUN } from '../../modules/abritamr/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow ABRITAMR {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    ABRITAMR_RUN(fasta)
    CSVTK_CONCAT(gather(ABRITAMR_RUN.out.summary, 'abritamr'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    summary_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.summary
    merged_summary_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    matches_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.matches
    partials_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.partials
    virulence_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.virulence
    amrfinder_tsv: Channel<Tuple<Map, Path>> = ABRITAMR_RUN.out.amrfinder

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.summary,
        ABRITAMR_RUN.out.matches,
        ABRITAMR_RUN.out.partials,
        ABRITAMR_RUN.out.virulence,
        ABRITAMR_RUN.out.amrfinder,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ABRITAMR_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
