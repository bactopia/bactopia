//
// ismapper - Identify insertion sites positions in bacterial genomes
//
nextflow.preview.types = true

include { ISMAPPER as ISMAPPER_MODULE } from '../../modules/ismapper/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow ISMAPPER {
    take:
    ch_reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]
    ch_reference: Channel<Tuple<Map, Path>> // channel: reference genome file
    ch_insertions: Channel<Tuple<Map, Path>> // channel: insertion sequences file

    main:
    ISMAPPER_MODULE(ch_reads, ch_reference, ch_insertions)

    emit:
    results: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.supplemental
    logs: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = ISMAPPER_MODULE.out.versions
}
