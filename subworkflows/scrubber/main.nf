//
// scrubber - Scrub human reads from FASTQ files
//
nextflow.preview.types = true

include { SRAHUMANSCRUBBER } from '../srahumanscrubber/main'
include { K2SCRUBBER       } from '../k2scrubber/main'
include { CSVTK_CONCAT     } from '../../modules/csvtk/concat/main'
include { flattenPaths     } from 'plugin/nf-bactopia'
include { gather           } from 'plugin/nf-bactopia'

workflow SCRUBBER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    use_srascrubber: Boolean

    main:
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>
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

    CSVTK_CONCAT(gather(ch_scrub_report, 'scrubber'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = ch_scrub_report
    special_tsv: Channel<Tuple<Map, Path>> = ch_special_report
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    scrubbed: Channel<Tuple<Map, Path>> = ch_scrubbed
    scrubbed_extra: Channel<Tuple<Map, Path>> = ch_scrubbed_extra

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_scrub_report,
        ch_scrubbed,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_versions,
        CSVTK_CONCAT.out.versions
    ])
}
