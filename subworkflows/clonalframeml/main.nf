//
// clonalframeml - Predict recomination events in bacterial genomes
//
include { IQTREE } from '../../modules/iqtree/main'
include { CLONALFRAMEML as CLONALFRAME } from '../../modules/clonalframeml/main'

workflow CLONALFRAMEML {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // Create a quick start tree
    IQTREE(alignment)
    ch_versions = ch_versions.mix(IQTREE.out.versions)
    ch_logs = ch_logs.mix(IQTREE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)

    // Run ClonalFrameML
    CLONALFRAME(IQTREE.out.aln_tree)
    ch_versions = ch_versions.mix(CLONALFRAME.out.versions)
    ch_logs = ch_logs.mix(CLONALFRAME.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CLONALFRAME.out.nf_logs)

    emit:
    masked_aln = CLONALFRAME.out.masked_aln
    logs = ch_logs // channel: [ val(meta), [ logs ] ]
    nf_logs = ch_nf_logs // channel: [ val(meta), [ nf_logs ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
