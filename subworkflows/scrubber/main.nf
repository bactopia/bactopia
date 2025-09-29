//
// scrubber - Scrub human reads from FASTQ files
//
include { SRAHUMANSCRUBBER } from '../srahumanscrubber/main'
include { K2SCRUBBER } from '../k2scrubber/main'

workflow SCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_scrub_report = Channel.empty()
    ch_scrubbed = Channel.empty()
    ch_scrubbed_extra = Channel.empty()

    if (params.use_srascrubber) {
        SRAHUMANSCRUBBER(reads)
        ch_versions = ch_versions.mix(SRAHUMANSCRUBBER.out.versions)
        ch_logs = ch_logs.mix(SRAHUMANSCRUBBER.out.logs)
        ch_scrub_report = ch_scrub_report.mix(SRAHUMANSCRUBBER.out.scrub_report)
        ch_scrubbed = ch_scrubbed.mix(SRAHUMANSCRUBBER.out.scrubbed)
        ch_scrubbed_extra = ch_scrubbed_extra.mix(SRAHUMANSCRUBBER.out.scrubbed_extra)
    } else {
        K2SCRUBBER(reads)
        ch_versions = ch_versions.mix(K2SCRUBBER.out.versions)
        ch_logs = ch_logs.mix(K2SCRUBBER.out.logs)
        ch_scrub_report = ch_scrub_report.mix(K2SCRUBBER.out.scrub_report)
        ch_scrubbed = ch_scrubbed.mix(K2SCRUBBER.out.scrubbed)
        ch_scrubbed_extra = ch_scrubbed_extra.mix(K2SCRUBBER.out.scrubbed_extra)
    }

    // Collect all nf_logs
    ch_nf_logs = Channel.empty()
    if (params.use_srascrubber) {
        ch_nf_logs = ch_nf_logs.mix(
            SRAHUMANSCRUBBER.out.nf_logs
        )
    } else {
        ch_nf_logs = ch_nf_logs.mix(
            K2SCRUBBER.out.nf_logs,
        )
    }

    emit:
    // Individual outputs
    tsv = ch_scrub_report
    scrubbed = ch_scrubbed
    scrubbed_extra = ch_scrubbed_extra

    // Generic aggregate outputs
    results = ch_scrub_report.mix(
        ch_scrubbed,
        ch_scrubbed_extra
    )
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions
}
