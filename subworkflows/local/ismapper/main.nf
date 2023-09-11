//
// ismapper - Identify insertion sites positions in bacterial genomes
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'ismapper')
options.args = [
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
options.subdir = params.run_name
options.logs_use_prefix = true
REFERENCE = params.reference ? file(params.reference) : []
INSERTIONS = params.insertions ? file(params.insertions) : []
include { ISMAPPER as ISMAPPER_MODULE } from '../../../modules/nf-core/ismapper/main' addParams( options: options )

workflow ISMAPPER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    ISMAPPER_MODULE(reads, REFERENCE, INSERTIONS)
    ch_versions = ch_versions.mix(ISMAPPER_MODULE.out.versions.first())

    emit:
    results = ISMAPPER_MODULE.out.results
    versions = ch_versions // channel: [ versions.yml ]
}
