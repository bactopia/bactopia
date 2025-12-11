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
 * @modules iqtree as iqtree_module
 *
 * @input aln
 * Channel containing aln data
 *
 * @output phylogeny Phylogeny
 * @output alignment Alignment
 * @output aln_tree  Aln Tree
 * @output results   Aggregated results channel containing all output files
 * @output logs      Aggregated logs channel containing all execution logs
 * @output nf_logs   Aggregated Nextflow execution logs from all processes
 * @output versions  Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { IQTREE as IQTREE_MODULE } from '../../modules/iqtree/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow IQTREE {
    take:
    aln: Channel<Tuple<Map, Set<Path>>>

    main:
    IQTREE_MODULE(aln)

    emit:
    // Individual outputs
    phylogeny: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.phylogeny
    alignment: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.alignment
    aln_tree: Channel<Tuple<Map, Path, Path>> = IQTREE_MODULE.out.aln_tree

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        IQTREE_MODULE.out.phylogeny,
        IQTREE_MODULE.out.alignment,
        IQTREE_MODULE.out.supplemental
    ])
    logs: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.versions
}
