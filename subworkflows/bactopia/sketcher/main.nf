//
// sketcher - Assortment of tools for sketching sequences
//
nextflow.preview.types = true

include { SKETCHER as SKETCHER_MODULE } from '../../../modules/bactopia/sketcher/main'

workflow SKETCHER {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    mash_db // channel: [ mash_db ]
    sourmash_db // channel: [ sourmash_db ]

    main:
    SKETCHER_MODULE(reads, mash_db, sourmash_db)

    emit:
    // Individual outputs
    sig = SKETCHER_MODULE.out.sig
    msh = SKETCHER_MODULE.out.msh
    mash = SKETCHER_MODULE.out.mash
    sourmash = SKETCHER_MODULE.out.sourmash

    // Generic aggregate outputs
    results = SKETCHER_MODULE.out.sig.mix(
        SKETCHER_MODULE.out.msh,
        SKETCHER_MODULE.out.mash,
        SKETCHER_MODULE.out.sourmash
    )
    logs = SKETCHER_MODULE.out.logs
    nf_logs = SKETCHER_MODULE.out.nf_begin.mix(
        SKETCHER_MODULE.out.nf_err,
        SKETCHER_MODULE.out.nf_log,
        SKETCHER_MODULE.out.nf_out,
        SKETCHER_MODULE.out.nf_run,
        SKETCHER_MODULE.out.nf_sh,
        SKETCHER_MODULE.out.nf_trace
    )
    versions = SKETCHER_MODULE.out.versions
}
