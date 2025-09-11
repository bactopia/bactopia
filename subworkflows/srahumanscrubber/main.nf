//
// srahumanscrubber - Scrub human reads from FASTQ files using SRA Human Scrubber
//
include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB } from '../../modules/srahumanscrubber/scrub/main'

workflow SRAHUMANSCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)
    ch_versions = ch_versions.mix(SRAHUMANSCRUBBER_SCRUB.out.versions)

    emit:
    scrubbed = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    scrubbed_extra = SRAHUMANSCRUBBER_SCRUB.out.scrubbed_extra
    scrub_report = SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    logs = SRAHUMANSCRUBBER_SCRUB.out.logs
    nf_logs = SRAHUMANSCRUBBER_SCRUB.out.nf_begin.mix(
        SRAHUMANSCRUBBER_SCRUB.out.nf_err,
        SRAHUMANSCRUBBER_SCRUB.out.nf_log,
        SRAHUMANSCRUBBER_SCRUB.out.nf_out,
        SRAHUMANSCRUBBER_SCRUB.out.nf_run,
        SRAHUMANSCRUBBER_SCRUB.out.nf_sh,
        SRAHUMANSCRUBBER_SCRUB.out.nf_trace
    )
    versions = ch_versions // channel: [ versions.yml ]
}
