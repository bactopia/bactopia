//
// scrubber - Scrub human reads from FASTQ files
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'scrubber')
options.is_module = params.wf == 'scrubber' ? true : false
options.args = ""
options.ignore = [".db"]

/*
include { SRAHUMANSCRUBBER_INITDB } from '../../../modules/nf-core/srahumanscrubber/initdb/main' addParams( )

if (params.wf == 'teton') {
    include { SRAHUMANSCRUBBER_SCRUB_TETON as SRAHUMANSCRUBBER_SCRUB } from '../../../modules/nf-core/srahumanscrubber/scrub/main' addParams( options: options )
} else if (params.wf == 'cleanyerreads') {
    include { SRAHUMANSCRUBBER_SCRUB_MAIN as SRAHUMANSCRUBBER_SCRUB } from '../../../modules/nf-core/srahumanscrubber/scrub/main' addParams( options: options )
} else {
    include { SRAHUMANSCRUBBER_SCRUB } from '../../../modules/nf-core/srahumanscrubber/scrub/main' addParams( options: options )
}
*/

include { K2SCRUBBER } from '../k2scrubber/main' addParams( options: options )

workflow SCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    K2SCRUBBER(reads)

    //SRAHUMANSCRUBBER_INITDB()
    //SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)
    ch_versions = ch_versions.mix(K2SCRUBBER.out.versions.first())

    emit:
    scrubbed = K2SCRUBBER.out.scrubbed
    tsv = K2SCRUBBER.out.scrub_report
    versions = ch_versions // channel: [ versions.yml ]
}
