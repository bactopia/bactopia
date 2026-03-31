/**
 * Construct maximum likelihood phylogenetic trees from alignments.
 *
 * This subworkflow uses [IQ-TREE](https://github.com/Cibiv/IQ-TREE) to build
 * maximum likelihood phylogenetic trees from multiple sequence alignments. IQ-TREE
 * implements fast and effective stochastic algorithms for phylogenetic inference,
 * including automatic model selection via ModelFinder. It produces phylogenetic trees
 * with bootstrap support and various supplementary files for tree visualization
 * and analysis.
 *
 * @status stable
 * @keywords phylogeny, maximum likelihood, tree, bootstrap, model selection
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation iqtree
 *
 * @modules iqtree
 *
 * @input record(meta, aln)
 * - `meta`: Groovy Map containing sample information
 * - `aln`: Multiple sequence alignment in FASTA format
 *
 * @output sample_outputs
 * - `aln`: Input multiple sequence alignment (passed through)
 * - `nwk`: Maximum-likelihood phylogenetic tree in Newick format
 * - `supplemental`: Detailed report, distance matrix, and model parameters
 *
 * @output run_outputs
 */
nextflow.preview.types = true

include { IQTREE as IQTREE_MODULE } from '../../modules/iqtree/main'

workflow IQTREE {
    take:
    alignment: Channel<Record>

    main:
    IQTREE_MODULE(alignment)

    emit:
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = IQTREE_MODULE.out
}
