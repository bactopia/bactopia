//
// scrubber - Scrub human reads from FASTQ files
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'scrubber')
options.is_module = params.wf == 'scrubber' ? true : false
options.args = ""
options.ignore = [".db", "EMPTY_EXTRA"]

if (params.use_srascrubber) {
    include { SRAHUMANSCRUBBER as SCRUBBER_METHOD } from '../srahumanscrubber/main' addParams( options: options )
} else {
    include { K2SCRUBBER as SCRUBBER_METHOD } from '../k2scrubber/main' addParams( options: options )
}

workflow SCRUBBER {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    SCRUBBER_METHOD(reads)
    ch_versions = ch_versions.mix(SCRUBBER_METHOD.out.versions.first())

    emit:
    scrubbed = SCRUBBER_METHOD.out.scrubbed
    scrubbed_extra = SCRUBBER_METHOD.out.scrubbed_extra
    tsv = SCRUBBER_METHOD.out.scrub_report
    versions = ch_versions // channel: [ versions.yml ]
}
