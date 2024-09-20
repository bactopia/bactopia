//
// teton - Taxonomic classification and estimated species abundances
//
include { SCRUBBER } from '../scrubber/main'
include { BRACKEN } from '../bracken/main'
include { CSVTK_JOIN } from '../../../modules/nf-core/csvtk/join/main' addParams( options: [publish_to_base: true, logs_subdir: 'amrfinderplus-proteins'] )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: 'amrfinderplus-proteins'] )

workflow TETON {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    if (params.skip_scrubber) {
        // Taxon Classification & Abundance
        BRACKEN(reads)
        ch_versions = ch_versions.mix(BRACKEN.out.versions)
    } else {
        // Remove host reads
        SCRUBBER(reads)
        ch_versions = ch_versions.mix(SCRUBBER.out.versions)

        // Taxon Classification & Abundance
        BRACKEN(SCRUBBER.out.scrubbed)
        ch_versions = ch_versions.mix(BRACKEN.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
