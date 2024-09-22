//
// srahumanscrubber - Scrub human reads from FASTQ files using SRA Human Scrubber
//

include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'srahumanscrubber')
options.is_module = params.wf == 'scrubber' ? true : false
options.args = ""

include { SRAHUMANSCRUBBER_INITDB } from '../../../modules/nf-core/srahumanscrubber/initdb/main' addParams( )
include { SRAHUMANSCRUBBER_SCRUB } from '../../../modules/nf-core/srahumanscrubber/scrub/main' addParams( options: options )

workflow SRAHUMANSCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)
    ch_versions = ch_versions.mix(SRAHUMANSCRUBBER_SCRUB.out.versions.first())

    emit:
    scrubbed = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    scrubbed_extra = SRAHUMANSCRUBBER_SCRUB.out.scrubbed_extra
    scrub_report = SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    versions = ch_versions // channel: [ versions.yml ]
}
