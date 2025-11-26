//
// srahumanscrubber - Scrub human reads from FASTQ files using SRA Human Scrubber
//
nextflow.preview.types = true

include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB  } from '../../modules/srahumanscrubber/scrub/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SRAHUMANSCRUBBER {
    take:
    reads: Channel<Tuple<Map, Path>> // channel: [ val(meta), [ reads ] ]

    main:
    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)

    emit:
    // Individual outputs
    scrubbed: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    scrubbed_extra: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrubbed_extra
    scrub_report: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    scrub_special_report: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrub_special_report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SRAHUMANSCRUBBER_SCRUB.out.scrubbed,
        SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    ])
    logs: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.versions
}
