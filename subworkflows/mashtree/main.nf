/**
 * Create phylogenetic trees using Mash distances.
 *
 * This subworkflow uses [Mashtree](https://github.com/lskatz/mashtree) to rapidly compare
 * whole genome sequence files and generate phylogenetic trees. It creates Mash sketches
 * of input genomes, calculates pairwise distances, and constructs a tree based on
 * the distance matrix.
 *
 * @status stable
 * @keywords phylogeny, tree, mash, distance, comparison
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation mashtree
 *
 * @modules mashtree
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output matrix          Pairwise distance matrix containing Mash distances between all samples
 * @output sketches        Mash sketch files for each input assembly
 * @output tree            Newick-format phylogenetic tree constructed from Mash distances
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow MASHTREE {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    MASHTREE_MODULE(gather(assembly, 'mashtree', 'fna'))

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
    logs: Channel<Tuple<Map, Path>> = flattenPaths([MASHTREE_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([MASHTREE_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([MASHTREE_MODULE.out.versions])
}
