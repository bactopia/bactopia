//
// gubbins - Rapid phylogenetic analysis of recombinant bacterial sequences
//
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
