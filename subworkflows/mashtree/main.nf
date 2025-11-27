//
// mashtree - Quickly create a tree using Mash distances
//
nextflow.preview.types = true

include { MASHTREE as MASHTREE_MODULE } from '../../modules/mashtree/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow MASHTREE {
    take:
    fasta: Channel<Tuple<Map, Set<Path>>>

    main:
    MASHTREE_MODULE(gather(fasta, 'mashtree', 'fna'))

    emit:
    // Individual outputs
    matrix: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.matrix
    sketches: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.sketches
    tree: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.tree

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MASHTREE_MODULE.out.matrix,
        MASHTREE_MODULE.out.sketches,
        MASHTREE_MODULE.out.tree
    ])
    logs: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = MASHTREE_MODULE.out.versions
}
