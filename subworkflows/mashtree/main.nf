//
// mashtree - Quickly create a tree using Mash distances
//
include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'

workflow MASHTREE {
    take:
    fasta // channel: [ [meta], [assemblies] ]

    main:
    fasta.collect{_meta, fna -> fna}.map{ fna -> [[id: 'mashtree'], fna]}.set{ ch_merge_fna }
    MASHTREE_MODULE(ch_merge_fna)

    emit:
    matrix = MASHTREE_MODULE.out.matrix
    sketches = MASHTREE_MODULE.out.sketches
    tree = MASHTREE_MODULE.out.tree
    logs = MASHTREE_MODULE.out.logs
    nf_logs = MASHTREE_MODULE.out.nf_begin.mix(
        MASHTREE_MODULE.out.nf_err,
        MASHTREE_MODULE.out.nf_log,
        MASHTREE_MODULE.out.nf_out,
        MASHTREE_MODULE.out.nf_run,
        MASHTREE_MODULE.out.nf_sh,
        MASHTREE_MODULE.out.nf_trace
    )
    versions = MASHTREE_MODULE.out.versions
}
