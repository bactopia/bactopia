/**
 * Detect and mask recombination events in bacterial phylogenies.
 *
 * This subworkflow uses [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML) to
 * detect and mask recombination events in bacterial phylogenies. It first builds a quick
 * phylogenetic tree using [IQ-TREE](https://github.com/iqtree/iqtree2), then identifies
 * recombination regions and creates a recombination-masked alignment. Finally, it
 * calculates SNP distances from the masked alignment using [snp-dists](https://github.com/tseemann/snp-dists).
 *
 * @status stable
 * @keywords recombination, phylogeny, masking, clonalframe, bacterial evolution
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation clonalframeml, iqtree, snpdists
 *
 * @subworkflows iqtree, snpdists
 * @modules clonalframeml
 *
 * @input record(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Core-genome alignment in FASTA format
 *
 * @output sample_outputs
 * - `masked_aln`: Recombination-masked alignment with detected recombination regions removed
 * - `emsim`: Uncertainty estimation results
 * - `em`: Final parameter estimates from the EM algorithm
 * - `status`: Predicted recombination events (importations)
 * - `newick`: Input tree with internal nodes labelled
 * - `fasta`: Reconstructed ancestral sequences
 * - `pos_ref`: Position cross-reference table
 *
 * @output iqtree_outputs
 * - `msa`: Input multiple sequence alignment (passed through)
 * - `phylogeny`: Quick-start maximum-likelihood phylogenetic tree
 *
 * @output snpdists_outputs
 * - `tsv`: Pairwise SNP distances from masked alignment in TSV format
 */
nextflow.preview.types = true

include { CLONALFRAMEML as CLONALFRAMEML_MODULE } from '../../modules/clonalframeml/main'
include { IQTREE                                } from '../iqtree/main'
include { SNPDISTS                              } from '../snpdists/main'
include { gather                                } from 'plugin/nf-bactopia'

workflow CLONALFRAMEML {
    take:
    alignment: Channel<Record>

    main:
    // Create a quick start tree
    IQTREE(gather(alignment, 'iqtree-fast', args: 'name: "iqtree-fast", process_name: "iqtree-fast"'))

    // Run ClonalFrameML - gather alignment and tree together
    ch_aln_tree = IQTREE.out.sample_outputs.collect{ r -> [r.msa, r.phylogeny] }.map{ msa, phylogeny ->
        [[name: "core-genome", process_name: "clonalframeml"], msa, phylogeny]
    }
    CLONALFRAMEML_MODULE(ch_aln_tree)

    // Per-sample SNP distances
    SNPDISTS(gather(CLONALFRAMEML_MODULE.out, 'core-genome.masked.distance', field: 'masked_aln', args: 'name: "core-genome.masked.distance", process_name: "snpdists-masked"'))

    emit:
    sample_outputs = CLONALFRAMEML_MODULE.out
    iqtree_outputs = IQTREE.out.sample_outputs
    snpdists_outputs = SNPDISTS.out.sample_outputs
}
