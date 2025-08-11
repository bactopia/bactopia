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
    ch_nf_logs = Channel.empty()

    GUBBINS_MODULE(alignment)
    ch_versions = ch_versions.mix(GUBBINS_MODULE.out.versions)
    ch_logs = ch_logs.mix(GUBBINS_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GUBBINS_MODULE.out.nf_logs)

    emit:
    masked_aln = GUBBINS_MODULE.out.masked_aln
    fasta = GUBBINS_MODULE.out.fasta
    gff = GUBBINS_MODULE.out.gff
    vcf = GUBBINS_MODULE.out.vcf
    phylip = GUBBINS_MODULE.out.phylip
    embl_predicted = GUBBINS_MODULE.out.embl_predicted
    embl_branch = GUBBINS_MODULE.out.embl_branch
    tree = GUBBINS_MODULE.out.tree
    tree_labelled = GUBBINS_MODULE.out.tree_labelled
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
