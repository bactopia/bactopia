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
    ch_nf_logs = Channel.empty()
    ch_scrubbed = Channel.empty()
    ch_scrubbed_extra = Channel.empty()
    ch_scrub_report = Channel.empty()

    if (params.use_srascrubber) {
        SRAHUMANSCRUBBER(reads)
        ch_versions = ch_versions.mix(SRAHUMANSCRUBBER.out.versions.first())
        ch_logs = ch_logs.mix(SRAHUMANSCRUBBER.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SRAHUMANSCRUBBER.out.nf_logs)
        ch_scrubbed = ch_scrubbed.mix(SRAHUMANSCRUBBER.out.scrubbed)
        ch_scrubbed_extra = ch_scrubbed_extra.mix(SRAHUMANSCRUBBER.out.scrubbed_extra)
        ch_scrub_report = ch_scrub_report.mix(SRAHUMANSCRUBBER.out.scrub_report)
    } else {
        K2SCRUBBER(reads)
        ch_versions = ch_versions.mix(K2SCRUBBER.out.versions.first())
        ch_logs = ch_logs.mix(K2SCRUBBER.out.logs)
        ch_nf_logs = ch_nf_logs.mix(K2SCRUBBER.out.nf_logs)
        ch_scrubbed = ch_scrubbed.mix(K2SCRUBBER.out.scrubbed)
        ch_scrubbed_extra = ch_scrubbed_extra.mix(K2SCRUBBER.out.scrubbed_extra)
        ch_scrub_report = ch_scrub_report.mix(K2SCRUBBER.out.scrub_report)
    }

    emit:
    scrubbed = ch_scrubbed
    scrubbed_extra = ch_scrubbed_extra
    tsv = ch_scrub_report
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
