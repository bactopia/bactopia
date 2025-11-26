//
// sketcher - Assortment of tools for sketching sequences
//
nextflow.preview.types = true

include { SKETCHER as SKETCHER_MODULE } from '../../../modules/bactopia/sketcher/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SKETCHER {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]
    mash_db: Channel<Tuple<Map, Path>> // channel: [ mash_db ]
    sourmash_db: Channel<Tuple<Map, Path>> // channel: [ sourmash_db ]

    main:
    SKETCHER_MODULE(reads, mash_db, sourmash_db)

    emit:
    // Individual outputs
    sig: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.sig
    msh: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.msh
    mash: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.mash
    sourmash: Channel<Tuple<Map, Path>> = SKETCHER_MODULE.out.sourmash

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SKETCHER_MODULE.out.sig,
        SKETCHER_MODULE.out.msh,
        SKETCHER_MODULE.out.mash,
        SKETCHER_MODULE.out.sourmash
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SKETCHER_MODULE.out.versions])
}
