//
// gubbins - Rapid phylogenetic analysis of recombinant bacterial sequences
//
include { GUBBINS as GUBBINS_MODULE } from '../../modules/gubbins/main'

workflow GUBBINS {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    GUBBINS_MODULE(alignment)
    ch_versions = ch_versions.mix(GUBBINS_MODULE.out.versions)
    ch_logs = ch_logs.mix(GUBBINS_MODULE.out.logs)

    emit:
    vcf = GUBBINS_MODULE.out.vcf
    embl_branch = GUBBINS_MODULE.out.embl_branch
    embl_predicted = GUBBINS_MODULE.out.embl_predicted
    fasta = GUBBINS_MODULE.out.fasta
    gff = GUBBINS_MODULE.out.gff
    masked_aln = GUBBINS_MODULE.out.masked_aln
    phylip = GUBBINS_MODULE.out.phylip
    tree = GUBBINS_MODULE.out.tree
    tree_labelled = GUBBINS_MODULE.out.tree_labelled
    logs = ch_logs
    nf_logs = GUBBINS_MODULE.out.nf_begin.mix(
        GUBBINS_MODULE.out.nf_err,
        GUBBINS_MODULE.out.nf_log,
        GUBBINS_MODULE.out.nf_out,
        GUBBINS_MODULE.out.nf_run,
        GUBBINS_MODULE.out.nf_sh,
        GUBBINS_MODULE.out.nf_trace
    )
    versions = ch_versions
}
