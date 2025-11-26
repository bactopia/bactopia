//
// iqtree - Phylogeny from a multiple sequence alignment using the maxium likelihood algorithm
//
nextflow.preview.types = true

include { IQTREE as IQTREE_MODULE } from '../../modules/iqtree/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow IQTREE {
    take:
    aln: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]

    main:
    IQTREE_MODULE(aln)

    emit:
    // Individual outputs
    phylogeny: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.phylogeny
    alignment: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.alignment
    aln_tree: Channel<Tuple<Map, Path, Path>> = IQTREE_MODULE.out.aln_tree

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        IQTREE_MODULE.out.phylogeny,
        IQTREE_MODULE.out.alignment,
        IQTREE_MODULE.out.supplemental
    ])
    logs: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = IQTREE_MODULE.out.versions
}
