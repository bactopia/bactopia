//
// scrubber - Scrub human reads from FASTQ files
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'srahumanscrubber')
options.is_module = params.wf == 'scrubber' ? true : false
options.args = ""
options.ignore = [".db"]

include { SRAHUMANSCRUBBER_INITDB } from '../../../modules/nf-core/srahumanscrubber/initdb/main' addParams( )

if (params.is_subworkflow) {
    include { SRAHUMANSCRUBBER_SCRUB } from '../../../modules/nf-core/srahumanscrubber/scrub/main' addParams( options: options )
} else {
    include { SRAHUMANSCRUBBER_SCRUB_MAIN as SRAHUMANSCRUBBER_SCRUB } from '../../../modules/nf-core/srahumanscrubber/scrub/main' addParams( options: options )
}


workflow SCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    if (params.download_scrubber) {
        SRAHUMANSCRUBBER_INITDB()
        ch_versions = ch_versions.mix(SRAHUMANSCRUBBER_INITDB.out.versions.first())

        SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)
    } else {
        SRAHUMANSCRUBBER_SCRUB(reads, file(params.scrubber_db))
    }

    ch_versions = ch_versions.mix(SRAHUMANSCRUBBER_SCRUB.out.versions.first())

    emit:
    scrubbed = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    versions = ch_versions // channel: [ versions.yml ]
}
