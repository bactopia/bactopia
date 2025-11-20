//
// scrubber - Scrub human reads from FASTQ files
//
include { SRAHUMANSCRUBBER } from '../srahumanscrubber/main'
include { K2SCRUBBER       } from '../k2scrubber/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'

workflow SCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    use_srascrubber

    main:
    ch_results = channel.empty()
    ch_logs = channel.empty()
    ch_nf_logs = channel.empty()
    ch_versions = channel.empty()
    ch_scrub_report = channel.empty()
    ch_special_report = channel.empty()
    ch_scrubbed = channel.empty()
    ch_scrubbed_extra = channel.empty()

    if (use_srascrubber) {
        SRAHUMANSCRUBBER(reads)
        ch_results = ch_results.mix(SRAHUMANSCRUBBER.out.results)
        ch_logs = ch_logs.mix(SRAHUMANSCRUBBER.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SRAHUMANSCRUBBER.out.nf_logs)
        ch_versions = ch_versions.mix(SRAHUMANSCRUBBER.out.versions)
        ch_scrub_report = ch_scrub_report.mix(SRAHUMANSCRUBBER.out.scrub_report)
        ch_special_report = ch_special_report.mix(SRAHUMANSCRUBBER.out.scrub_special_report)
        ch_scrubbed = ch_scrubbed.mix(SRAHUMANSCRUBBER.out.scrubbed)
        ch_scrubbed_extra = ch_scrubbed_extra.mix(SRAHUMANSCRUBBER.out.scrubbed_extra)
    } else {
        K2SCRUBBER(reads)
        ch_results = ch_results.mix(K2SCRUBBER.out.results)
        ch_logs = ch_logs.mix(K2SCRUBBER.out.logs)
        ch_nf_logs = ch_nf_logs.mix(K2SCRUBBER.out.nf_logs)
        ch_versions = ch_versions.mix(K2SCRUBBER.out.versions)
        ch_scrub_report = ch_scrub_report.mix(K2SCRUBBER.out.scrub_report)
        ch_special_report = ch_special_report.mix(K2SCRUBBER.out.scrub_special_report)
        ch_scrubbed = ch_scrubbed.mix(K2SCRUBBER.out.scrubbed)
        ch_scrubbed_extra = ch_scrubbed_extra.mix(K2SCRUBBER.out.scrubbed_extra)
    }

    // Merge results
    ch_scrub_report.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'scrubber'], tsv]}.set{ ch_merge_sccmec }
    CSVTK_CONCAT(ch_merge_sccmec, 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv = ch_scrub_report
    special_tsv = ch_special_report
    merged_tsv = CSVTK_CONCAT.out.csv
    scrubbed = ch_scrubbed
    scrubbed_extra = ch_scrubbed_extra

    // Generic aggregate outputs
    results = ch_scrub_report.mix(
        ch_scrubbed,
        CSVTK_CONCAT.out.csv
    )
    logs = ch_logs.mix(
        CSVTK_CONCAT.out.logs
    )
    nf_logs = ch_nf_logs.mix(
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace
    )
    versions = ch_versions.mix(
        CSVTK_CONCAT.out.versions
    )
}
