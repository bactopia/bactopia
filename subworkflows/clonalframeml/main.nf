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


workflow CLONALFRAMEML {
    take:
    alignment: Channel<Record>

    main:
    // Create a quick start tree
    IQTREE(alignment.map { r ->
        record(_meta: [name: r._meta.name, process_name: 'iqtree-fast'], msa: r.alignment)
    })

    // Run ClonalFrameML
    CLONALFRAMEML_MODULE(IQTREE.out.sample_outputs.map { r ->
        record(_meta: [name: "core-genome", process_name: "clonalframeml"], msa: r.msa, newick: r.phylogeny)
    })

    // Per-sample SNP distances
    SNPDISTS(CLONALFRAMEML_MODULE.out.map { r ->
        record(_meta: [name: 'core-genome.masked.distance', process_name: 'snpdists-masked'], msa: r.masked_aln)
    })

    emit:
    // Published outputs
    sample_outputs = CLONALFRAMEML_MODULE.out
    iqtree_outputs = IQTREE.out.sample_outputs
    snpdists_outputs = SNPDISTS.out.sample_outputs
}
