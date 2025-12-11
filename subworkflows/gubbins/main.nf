/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @subworkflows snpdists
 * @modules gubbins as gubbins_module
 *
 * @input alignment
 * Channel containing alignment data
 *
 * @output masked_aln     Masked Aln
 * @output fasta          Fasta
 * @output gff            Gff
 * @output vcf            Vcf
 * @output stats          Stats
 * @output phylip         Phylip
 * @output embl_predicted Embl Predicted
 * @output embl_branch    Embl Branch
 * @output tree           Tree
 * @output tree_labelled  Tree Labelled
 * @output bootstrap_tree Bootstrap Tree
 * @output tsv            Tsv
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution logs from all processes
 * @output versions       Aggregated version information from all executed tools
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
    masked_aln: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.masked_aln
    fasta: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.fasta
    gff: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.gff
    vcf: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.vcf
    stats: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.stats
    phylip: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.phylip
    embl_predicted: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.embl_predicted
    embl_branch: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.embl_branch
    tree: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.tree
    tree_labelled: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.tree_labelled
    bootstrap_tree: Channel<Tuple<Map, Path>> = GUBBINS_MODULE.out.bootstrap_tree
    tsv: Channel<Tuple<Map, Path>> = SNPDISTS.out.tsv

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
