//
// ismapper - Identify insertion sites positions in bacterial genomes
//
ismapper_args = [
    params.ismap_all ? "-all" : "",
    "--min_clip ${params.min_clip}",
    "--max_clip ${params.max_clip}",
    "--cutoff ${params.cutoff}",
    "--novel_gap_size ${params.novel_gap_size}",
    "--min_range ${params.min_range}",
    "--max_range ${params.max_range}",
    "--merging ${params.merging}",
    "--T ${params.ismap_minqual}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { ISMAPPER as ISMAPPER_MODULE } from '../../../modules/nf-core/modules/ismapper/main' addParams( options: [ args: "${ismapper_args}", is_module: true] )

workflow ISMAPPER {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    reference // file
    query // file

    main:
    ch_versions = Channel.empty()

    ISMAPPER_MODULE(reads, reference, query)
    ch_versions = ch_versions.mix(ISMAPPER_MODULE.out.versions.first())

    emit:
    results = ISMAPPER_MODULE.out.results
    versions = ch_versions // channel: [ versions.yml ]
}
