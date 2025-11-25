//
// iqtree - Phylogeny from a multiple sequence alignment using the maxium likelihood algorithm
//
nextflow.preview.types = true

include { IQTREE as IQTREE_MODULE } from '../../modules/iqtree/main'

workflow IQTREE {
    take:
    aln // channel: [ val(meta), [ reads ] ]

    main:
    IQTREE_MODULE(aln)

    emit:
    // Individual outputs
    phylogeny = IQTREE_MODULE.out.phylogeny
    alignment = IQTREE_MODULE.out.alignment
    aln_tree = IQTREE_MODULE.out.aln_tree

    // Generic aggregate outputs
    results = IQTREE_MODULE.out.phylogeny.mix(
        IQTREE_MODULE.out.alignment,
        IQTREE_MODULE.out.supplemental
    )
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
