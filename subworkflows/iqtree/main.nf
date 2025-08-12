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
    IQTREE_MODULE(aln)
    ch_versions = ch_versions.mix(IQTREE_MODULE.out.versions.first())
    ch_logs = ch_logs.mix(IQTREE_MODULE.out.logs)

    emit:
    phylogeny = IQTREE_MODULE.out.phylogeny
    logs = ch_logs
    nf_logs = IQTREE_MODULE.out.nf_begin.mix(
        IQTREE_MODULE.out.nf_err,
        IQTREE_MODULE.out.nf_log,
        IQTREE_MODULE.out.nf_out,
        IQTREE_MODULE.out.nf_run,
        IQTREE_MODULE.out.nf_sh,
        IQTREE_MODULE.out.nf_trace
    )
    versions = ch_versions
}
