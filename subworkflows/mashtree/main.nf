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
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tree`: Phylogenetic tree in Newick format
 * - `tsv`: Pairwise distance matrix
 * - `sketches`: Individual Mash sketch files
 */
nextflow.preview.types = true

include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'
include { gather                      } from 'plugin/nf-bactopia'

workflow MASHTREE {
    take:
    assembly: Channel<Record>

    main:
    MASHTREE_MODULE(gather(assembly, 'assembly', [name: 'mashtree']))

    emit:
    sample_outputs = MASHTREE_MODULE.out
}
