/**
 * Detect and filter recombination regions in bacterial alignments.
 *
 * This subworkflow uses [Gubbins](https://github.com/nickjcroucher/gubbins) (Globally
 * Optimised Bacterial Phylogenomic analysis) to identify recombination regions in
 * bacterial core-genome alignments. It iteratively filters out recombination to
 * produce a recombination-free phylogeny, then calculates SNP distances from the
 * masked alignment. Gubbins is essential for accurate phylogenetic reconstruction
 * of recombining bacterial species.
 *
 * @status stable
 * @keywords recombination, phylogeny, filter, snp, core-genome
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation
 * @citation gubbins, snpdists
 *
 * @subworkflows snpdists
 * @modules gubbins
 *
 * @input record(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA format
 *
 * @output sample_outputs
 * - `masked_aln`: Recombination-masked alignment in FASTA format
 * - `fasta`: Concatenated alignment before masking in FASTA format
 * - `gff`: GFF file containing recombination region coordinates
 * - `vcf`: VCF file containing SNPs filtered by Gubbins
 * - `stats`: Summary statistics of the Gubbins analysis
 * - `phylip`: Recombination-masked alignment in PHYLIP format
 * - `embl_predicted`: Recombination predictions in EMBL format
 * - `embl_branch`: Branch-specific recombination in EMBL format
 * - `tree`: Maximum likelihood tree from filtered SNPs in Newick format
 * - `tree_labelled`: Annotated tree with node labels in Newick format
 * - `bootstrap_tree`: Bootstrapped phylogenetic tree in Newick format
 *
 * @output snpdists_outputs
 * - `tsv`: Pairwise SNP distances from masked alignment in TSV format
 */
nextflow.preview.types = true

include { GUBBINS as GUBBINS_MODULE } from '../../modules/gubbins/main'
include { SNPDISTS                  } from '../snpdists/main'
include { gather                    } from 'plugin/nf-bactopia'

workflow GUBBINS {
    take:
    alignment: Channel<Record>

    main:
    GUBBINS_MODULE(alignment)
    SNPDISTS(gather(GUBBINS_MODULE.out, 'core-snp.masked.distance', field: 'masked_aln', args: 'name: "core-snp.masked.distance", process_name: "snpdists-masked"'))

    emit:
    sample_outputs = GUBBINS_MODULE.out
    snpdists_outputs = SNPDISTS.out.sample_outputs
}
