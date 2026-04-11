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
 * @modules clonalframeml
 * @subworkflows iqtree, snpdists
 *
 * @input record(meta, aln)
 * - `meta`: Groovy Record containing sample information
 * - `aln`: Core-genome alignment in FASTA format
 *
 * @output sample_outputs
 *
 * @output run_outputs
 * - `masked_aln`: Recombination-masked alignment with detected recombination regions removed
 * - `emsim`: Uncertainty estimation results
 * - `em`: Final parameter estimates from the EM algorithm
 * - `status`: Predicted recombination events (importations)
 * - `nwk`: Input tree with internal nodes labelled
 * - `fasta`: Reconstructed ancestral sequences
 * - `pos_ref`: Position cross-reference table
 * - `aln`: Input multiple sequence alignment (passed through from IQ-TREE)
 * - `nwk`: Quick-start maximum-likelihood phylogenetic tree from IQ-TREE
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
    ch_iqtree = IQTREE(alignment.map { r ->
        record(meta: record(name: r.meta.name, process_name: 'iqtree-fast'), aln: r.aln)
    })

    // Run ClonalFrameML
    ch_clonalframeml = CLONALFRAMEML_MODULE(ch_iqtree.run_outputs.map { r ->
        record(meta: record(name: "core-genome", process_name: "clonalframeml"), aln: r.aln, nwk: r.nwk)
    })

    // Per-sample SNP distances
    ch_snpdists = SNPDISTS(ch_clonalframeml.map { r ->
        record(meta: record(name: 'core-genome.masked.distance', process_name: 'snpdists-masked'), aln: r.masked_aln)
    })

    emit: // bactopia-lint: ignore S005, S010
    // Downstream inputs
    alignment = ch_clonalframeml.map { r ->
        record(meta: record(name: "core-genome", process_name: "iqtree"), aln: r.masked_aln)
    }
    // Published outputs
    sample_outputs = channel.empty()
    run_outputs = ch_clonalframeml.mix(ch_iqtree.run_outputs).mix(ch_snpdists.run_outputs)
}
