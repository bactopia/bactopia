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
 * @input tuple(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA format
 *
 * @output masked_aln          Recombination-masked alignment in FASTA format
 * @output fasta              Concatenated alignment before masking in FASTA format
 * @output gff                GFF file containing recombination region coordinates
 * @output vcf                VCF file containing SNPs filtered by Gubbins
 * @output stats              Summary statistics of the Gubbins analysis
 * @output phylip             Recombination-masked alignment in PHYLIP format
 * @output embl_predicted     Recombination predictions in EMBL format
 * @output embl_branch        Branch-specific recombination in EMBL format
 * @output tree               Maximum likelihood tree from filtered SNPs in Newick format
 * @output tree_labelled      Annotated tree with node labels in Newick format
 * @output bootstrap_tree     Bootstrapped phylogenetic tree in Newick format
 * @output tsv                Pairwise SNP distances from masked alignment in TSV format
 * @output results            Aggregated results channel containing all output files
 * @output logs               Aggregated logs channel containing all execution logs
 * @output nf_logs            Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions           Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { GUBBINS as GUBBINS_MODULE } from '../../modules/gubbins/main'
include { SNPDISTS                  } from '../snpdists/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow GUBBINS {
    take:
    alignment: Channel<Tuple<Map, Set<Path>>>

    main:
    GUBBINS_MODULE(alignment)
    SNPDISTS(gather(GUBBINS_MODULE.out.masked_aln, 'core-snp.masked.distance', 'name: "core-snp.masked.distance", process_name: "snpdists-masked"'))

    emit:
    // Individual outputs
    masked_aln: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.masked_aln
    fasta: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.fasta
    gff: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.gff
    vcf: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.vcf
    stats: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.stats
    phylip: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.phylip
    embl_predicted: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.embl_predicted
    embl_branch: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.embl_branch
    tree: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.tree
    tree_labelled: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.tree_labelled
    bootstrap_tree: Channel<Tuple<Map, Set<Path>>> = GUBBINS_MODULE.out.bootstrap_tree
    tsv: Channel<Tuple<Map, Set<Path>>> = SNPDISTS.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GUBBINS_MODULE.out.masked_aln,
        GUBBINS_MODULE.out.fasta,
        GUBBINS_MODULE.out.gff,
        GUBBINS_MODULE.out.vcf,
        GUBBINS_MODULE.out.stats,
        GUBBINS_MODULE.out.phylip,
        GUBBINS_MODULE.out.embl_predicted,
        GUBBINS_MODULE.out.embl_branch,
        GUBBINS_MODULE.out.tree,
        GUBBINS_MODULE.out.tree_labelled,
        GUBBINS_MODULE.out.bootstrap_tree,
        SNPDISTS.out.results
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GUBBINS_MODULE.out.logs,
        SNPDISTS.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GUBBINS_MODULE.out.nf_logs,
        SNPDISTS.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GUBBINS_MODULE.out.versions,
        SNPDISTS.out.versions
    ])
}
