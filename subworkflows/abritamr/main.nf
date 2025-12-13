/**
 * Identify antimicrobial resistance genes using AMRFinderPlus.
 *
 * This subworkflow uses [abriTAMR](https://github.com/MDU-PHL/abritamr) to identify
 * antimicrobial resistance genes in bacterial genomes. It runs AMRFinderPlus on
 * each sample and collates the results into functional classes, producing detailed
 * reports on resistance genes, partial matches, and virulence factors.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, amr, amrfinderplus, classification
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation abritamr
 *
 * @modules abritamr_run, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output summary_tsv        Per-sample TSV files with AMR gene summaries
 * @output merged_summary_tsv Consolidated TSV file containing AMR summaries from all samples
 * @output matches_tsv        Per-sample TSV files with complete AMR gene matches
 * @output partials_tsv       Per-sample TSV files with partial AMR gene matches
 * @output virulence_tsv      Per-sample TSV files with virulence gene predictions
 * @output amrfinder_tsv      Per-sample raw AMRFinderPlus output files
 * @output results            Aggregated results channel containing all output files
 * @output logs               Aggregated logs channel containing all execution logs
 * @output nf_logs            Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions           Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { ABRITAMR_RUN } from '../../modules/abritamr/run/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow ABRITAMR {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    ABRITAMR_RUN(assembly)
    CSVTK_CONCAT(gather(ABRITAMR_RUN.out.summary, 'abritamr'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    summary_tsv: Channel<Tuple<Map, Set<Path>>> = ABRITAMR_RUN.out.summary
    merged_summary_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    matches_tsv: Channel<Tuple<Map, Set<Path>>> = ABRITAMR_RUN.out.matches
    partials_tsv: Channel<Tuple<Map, Set<Path>>> = ABRITAMR_RUN.out.partials
    virulence_tsv: Channel<Tuple<Map, Set<Path>>> = ABRITAMR_RUN.out.virulence
    amrfinder_tsv: Channel<Tuple<Map, Set<Path>>> = ABRITAMR_RUN.out.amrfinder

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
