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
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Pre-gathered assembled contigs in FASTA format (multiple genomes)
 *
 * @output sample_outputs
 *
 * @output run_outputs
 * - `nwk`: Phylogenetic tree in Newick format
 * - `tsv`: Pairwise distance matrix
 * - `sketches`: Individual Mash sketch files
 */
nextflow.enable.types = true

include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'

workflow MASHTREE {
    take:
    assemblies: Channel<Record>

    main:
    ch_mashtree = MASHTREE_MODULE(assemblies)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_mashtree
}
