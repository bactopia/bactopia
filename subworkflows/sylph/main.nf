//
// sylph - Taxonomic profiling by abundance-corrected minhash
//
nextflow.preview.types = true

include { SYLPH_PROFILE } from '../../modules/sylph/profile/main'
include { flattenPaths  } from 'plugin/nf-bactopia'
include { gather        } from 'plugin/nf-bactopia'

workflow SYLPH {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]
    database: Channel<Tuple<Map, Path>>

    main:
    SYLPH_PROFILE(reads, database)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.tsv
    logs: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SYLPH_PROFILE.out.versions
}
