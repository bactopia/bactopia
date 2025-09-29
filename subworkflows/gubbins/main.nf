//
// gubbins - Rapid phylogenetic analysis of recombinant bacterial sequences
//
include { GUBBINS as GUBBINS_MODULE } from '../../modules/gubbins/main'

workflow GUBBINS {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    GUBBINS_MODULE(alignment)

    emit:
    // Individual outputs
    masked_aln = GUBBINS_MODULE.out.masked_aln
    fasta = GUBBINS_MODULE.out.fasta
    gff = GUBBINS_MODULE.out.gff
    vcf = GUBBINS_MODULE.out.vcf
    stats = GUBBINS_MODULE.out.stats
    phylip = GUBBINS_MODULE.out.phylip
    embl_predicted = GUBBINS_MODULE.out.embl_predicted
    embl_branch = GUBBINS_MODULE.out.embl_branch
    tree = GUBBINS_MODULE.out.tree
    tree_labelled = GUBBINS_MODULE.out.tree_labelled
    bootstrap_tree = GUBBINS_MODULE.out.bootstrap_tree

    // Generic aggregate outputs
    results = GUBBINS_MODULE.out.masked_aln.mix(
        GUBBINS_MODULE.out.fasta,
        GUBBINS_MODULE.out.gff,
        GUBBINS_MODULE.out.vcf,
        GUBBINS_MODULE.out.stats,
        GUBBINS_MODULE.out.phylip,
        GUBBINS_MODULE.out.embl_predicted,
        GUBBINS_MODULE.out.embl_branch,
        GUBBINS_MODULE.out.tree,
        GUBBINS_MODULE.out.tree_labelled,
        GUBBINS_MODULE.out.bootstrap_tree
    )
    logs = GUBBINS_MODULE.out.logs
    nf_logs = GUBBINS_MODULE.out.nf_begin.mix(
        GUBBINS_MODULE.out.nf_err,
        GUBBINS_MODULE.out.nf_log,
        GUBBINS_MODULE.out.nf_out,
        GUBBINS_MODULE.out.nf_run,
        GUBBINS_MODULE.out.nf_sh,
        GUBBINS_MODULE.out.nf_trace
    )
    versions = GUBBINS_MODULE.out.versions
}
