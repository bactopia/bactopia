//
// iqtree - Phylogeny from a multiple sequence alignment using the maxium likelihood algorithm
//
include { IQTREE as IQTREE_MODULE } from '../../modules/iqtree/main'

workflow IQTREE {
    take:
    aln // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    IQTREE_MODULE(aln)
    ch_versions = ch_versions.mix(IQTREE_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(IQTREE_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(IQTREE_MODULE.out.nf_logs)

    emit:
    phylogeny = IQTREE_MODULE.out.phylogeny
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
