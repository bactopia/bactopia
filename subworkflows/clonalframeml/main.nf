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
    // Create a quick start tree
    IQTREE(alignment)
    ch_versions = ch_versions.mix(IQTREE.out.versions)
    ch_logs = ch_logs.mix(IQTREE.out.logs)
    // Run ClonalFrameML
    CLONALFRAME(IQTREE.out.aln_tree)
    ch_versions = ch_versions.mix(CLONALFRAME.out.versions)
    ch_logs = ch_logs.mix(CLONALFRAME.out.logs)

    emit:

    emit:
    masked_aln = CLONALFRAME.out.masked_aln
    logs = ch_logs
    nf_logs = CLONALFRAME.out.nf_begin.mix(
        CLONALFRAME.out.nf_err,
        CLONALFRAME.out.nf_log,
        CLONALFRAME.out.nf_out,
        CLONALFRAME.out.nf_run,
        CLONALFRAME.out.nf_sh,
        CLONALFRAME.out.nf_trace,
        IQTREE.out.nf_begin,
        IQTREE.out.nf_err,
        IQTREE.out.nf_log,
        IQTREE.out.nf_out,
        IQTREE.out.nf_run,
        IQTREE.out.nf_sh,
        IQTREE.out.nf_trace
    )
    versions = ch_versions
}
