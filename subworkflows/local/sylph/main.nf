//
// sylph - Taxonomic profiling by abundance-corrected minhash
//
include { initOptions } from '../../../lib/nf/functions'
options  = initOptions(params.containsKey("options") ? params.options : [:], 'sylph')
options .args = [
    params.sylph_estimate_unknown ? "--estimate-unknown" : "",
    "-k ${params.sylph_k}",
    "-c ${params.sylph_subsample_rate}",
    "--min-spacing ${params.sylph_min_spacing}",
    "--minimum-ani ${params.sylph_min_ani}",
    "--min-number-kmers ${params.sylph_min_kmers}",
    "--min-count-correct ${params.sylph_min_correct}",
    params.sylph_opts ? "${params.sylph_opts}" : "",
].join(' ').replaceAll("\\s{2,}", " ").trim()
DATABASE = params.sylph_db ? file(params.sylph_db) : []

include { SYLPH_PROFILE }  from '../../../modules/nf-core/sylph/profile/main' addParams( options: options )

workflow SYLPH {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    // Run sylph profile
    SYLPH_PROFILE(reads, DATABASE)
    ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions.first())

    emit:
    tsv = SYLPH_PROFILE.out.tsv
    versions = ch_versions // channel: [ versions.yml ]
}
