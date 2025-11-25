//
// sylph - Taxonomic profiling by abundance-corrected minhash
//
nextflow.preview.types = true

include { SYLPH_PROFILE } from '../../modules/sylph/profile/main'

workflow SYLPH {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    database

    main:
    SYLPH_PROFILE(reads, database)

    emit:
    // Individual outputs
    tsv = SYLPH_PROFILE.out.tsv

    // Generic aggregate outputs
    results = SYLPH_PROFILE.out.tsv
    logs = SYLPH_PROFILE.out.logs
    nf_logs = SYLPH_PROFILE.out.nf_begin.mix(
        SYLPH_PROFILE.out.nf_err,
        SYLPH_PROFILE.out.nf_log,
        SYLPH_PROFILE.out.nf_out,
        SYLPH_PROFILE.out.nf_run,
        SYLPH_PROFILE.out.nf_sh,
        SYLPH_PROFILE.out.nf_trace
    )
    versions = SYLPH_PROFILE.out.versions
}
