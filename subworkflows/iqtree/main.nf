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
 * @input tuple(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA format
 *
 * @output phylogeny        Maximum likelihood tree in Newick format with bootstrap support
 * @output alignment        Processed alignment file used for tree construction
 * @output aln_tree         Combined alignment and tree file for visualization
 * @output results          Aggregated results channel containing all output files
 * @output logs             Aggregated logs channel containing all execution logs
 * @output nf_logs          Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions         Aggregated version information from all executed tools
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
    phylogeny: Channel<Tuple<Map, Set<Path>>> = IQTREE_MODULE.out.phylogeny
    alignment: Channel<Tuple<Map, Set<Path>>> = IQTREE_MODULE.out.alignment
    aln_tree: Channel<Tuple<Map, Set<Path>, Set<Path>>> = IQTREE_MODULE.out.aln_tree

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        IQTREE_MODULE.out.phylogeny,
        IQTREE_MODULE.out.alignment,
        IQTREE_MODULE.out.supplemental
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([IQTREE_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([IQTREE_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([IQTREE_MODULE.out.versions])
}
