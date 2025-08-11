//
// mashtree - Quickly create a tree using Mash distances
//
include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'

workflow MASHTREE {
    take:
    fasta // channel: [ val(meta), [ assemblies ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    MASHTREE_MODULE(fasta)
    ch_versions = ch_versions.mix(MASHTREE_MODULE.out.versions)
    ch_logs = ch_logs.mix(MASHTREE_MODULE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MASHTREE_MODULE.out.nf_logs)

    emit:
    tree = MASHTREE_MODULE.out.tree
    matrix = MASHTREE_MODULE.out.matrix
    sketches = MASHTREE_MODULE.out.sketches
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions
}
