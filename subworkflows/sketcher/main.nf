//
// sketcher - Assortment of tools for sketching sequences
//
include { SKETCHER as SKETCHER_MODULE } from '../../modules/bactopia/sketcher/main'

workflow SKETCHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    mash_db // channel: [ mash_db ]
    sourmash_db // channel: [ sourmash_db ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    // Sketch FASTQs
    SKETCHER_MODULE(reads, mash_db, sourmash_db)
    ch_versions = ch_versions.mix(SKETCHER_MODULE.out.versions)
    ch_logs = ch_logs.mix(SKETCHER_MODULE.out.logs)

    emit:

    emit:
    sig = SKETCHER_MODULE.out.sig
    logs = ch_logs
    nf_logs = SKETCHER_MODULE.out.nf_begin.mix(
        SKETCHER_MODULE.out.nf_err,
        SKETCHER_MODULE.out.nf_log,
        SKETCHER_MODULE.out.nf_out,
        SKETCHER_MODULE.out.nf_run,
        SKETCHER_MODULE.out.nf_sh,
        SKETCHER_MODULE.out.nf_trace
    )
    versions = ch_versions
}
