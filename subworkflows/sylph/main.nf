//
// sylph - Taxonomic profiling by abundance-corrected minhash
//
include { SYLPH_PROFILE } from '../../modules/sylph/profile/main'

workflow SYLPH {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    DATABASE = params.sylph_db ? file(params.sylph_db) : []

    // Run sylph profile
    SYLPH_PROFILE(reads, DATABASE)
    ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions.first()

    emit:
    tsv = SYLPH_PROFILE.out.tsv
    logs = SYLPH_PROFILE.out.logs
    nf_logs = SYLPH_PROFILE.out.nf_begin.mix(
        SYLPH_PROFILE.out.nf_err,
        SYLPH_PROFILE.out.nf_log,
        SYLPH_PROFILE.out.nf_out,
        SYLPH_PROFILE.out.nf_run,
        SYLPH_PROFILE.out.nf_sh,
        SYLPH_PROFILE.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
