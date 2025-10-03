//
// srahumanscrubber - Scrub human reads from FASTQ files using SRA Human Scrubber
//
include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB } from '../../modules/srahumanscrubber/scrub/main'

workflow SRAHUMANSCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)

    emit:
    // Individual outputs
    scrubbed = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    scrubbed_extra = SRAHUMANSCRUBBER_SCRUB.out.scrubbed_extra
    scrub_report = SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    scrub_special_report = SRAHUMANSCRUBBER_SCRUB.out.scrub_special_report

    // Generic aggregate outputs
    results = SRAHUMANSCRUBBER_SCRUB.out.scrubbed.mix(
        SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    )
    logs = SRAHUMANSCRUBBER_SCRUB.out.logs
    nf_logs = SRAHUMANSCRUBBER_SCRUB.out.nf_begin.mix(
        SRAHUMANSCRUBBER_SCRUB.out.nf_err,
        SRAHUMANSCRUBBER_SCRUB.out.nf_log,
        SRAHUMANSCRUBBER_SCRUB.out.nf_out,
        SRAHUMANSCRUBBER_SCRUB.out.nf_run,
        SRAHUMANSCRUBBER_SCRUB.out.nf_sh,
        SRAHUMANSCRUBBER_SCRUB.out.nf_trace
    )
    versions = SRAHUMANSCRUBBER_SCRUB.out.versions
}
