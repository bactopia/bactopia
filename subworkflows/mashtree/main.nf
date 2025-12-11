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
 * @modules mashtree as mashtree_module
 *
 * @input fasta
 * Channel containing fasta data
 *
 * @output matrix   Matrix
 * @output sketches Sketches
 * @output tree     Tree
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow MASHTREE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    MASHTREE_MODULE(gather(fasta, 'mashtree', 'fna'))

    emit:
    // Individual outputs
    matrix: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.matrix
    sketches: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.sketches
    tree: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.tree

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MASHTREE_MODULE.out.matrix,
        MASHTREE_MODULE.out.sketches,
        MASHTREE_MODULE.out.tree
    ])
    logs: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.versions
}
