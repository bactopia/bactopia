//
// SUBWORKFLOW_NAME - SUBWORKFLOW_DESCRIPTION
//

params.OPTS = [:]

include { MODULE } from '../../../path/to/module/main' addParams( options: params.OPTS )

workflow SUBWORKFLOW {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE DESCRIPTION
    //
    MODULE ( INPUTS )
    ch_versions = ch_versions.mix(MODULE.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
