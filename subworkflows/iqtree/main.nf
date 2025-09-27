//
// iqtree - Phylogeny from a multiple sequence alignment using the maxium likelihood algorithm
//
include { IQTREE as IQTREE_MODULE } from '../../modules/iqtree/main'

workflow IQTREE {
    take:
    aln // channel: [ val(meta), [ reads ] ]

    main:
    IQTREE_MODULE(aln)

    emit:
    phylogeny = IQTREE_MODULE.out.phylogeny
    aln_tree = IQTREE_MODULE.out.aln_tree
    results = IQTREE_MODULE.out.results.mix(IQTREE_MODULE.out.phylogeny)
    logs = IQTREE_MODULE.out.logs
    nf_logs = IQTREE_MODULE.out.nf_begin.mix(
        IQTREE_MODULE.out.nf_err,
        IQTREE_MODULE.out.nf_log,
        IQTREE_MODULE.out.nf_out,
        IQTREE_MODULE.out.nf_run,
        IQTREE_MODULE.out.nf_sh,
        IQTREE_MODULE.out.nf_trace
    )
    versions = IQTREE_MODULE.out.versions
}
