//
// pneumocat - Assign capsular type to Streptococcus pneumoniae from sequence reads
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'pneumocat')

include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../../modules/nf-core/pneumocat/main' addParams( options: options )

workflow PNEUMOCAT {
    take:
    fastq // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    PNEUMOCAT_MODULE(fastq)
    ch_versions = ch_versions.mix(PNEUMOCAT_MODULE.out.versions.first())

    emit:
    xml = PNEUMOCAT_MODULE.out.xml
    versions = ch_versions
}
