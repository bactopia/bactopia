//
// teton - Taxonomic classification and estimated species abundances
//
include { SCRUBBER } from '../scrubber/main'
include { BRACKEN } from '../bracken/main'
include { BACTOPIA_SAMPLESHEET } from '../../../modules/local/bactopia/teton/main'
include { CSVTK_JOIN } from '../../../modules/nf-core/csvtk/join/main' addParams( options: [publish_to_base: true, logs_subdir: 'amrfinderplus-proteins'] )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: 'amrfinderplus-proteins'] )

workflow TETON {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_reads = Channel.empty()

    // Remove host reads
    SCRUBBER(reads)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)

    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed)
    ch_versions = ch_versions.mix(BRACKEN.out.versions)

    // Determine genome size and create sample sheet
    BACTOPIA_SAMPLESHEET(BRACKEN.out.classification)
    ch_versions = ch_versions.mix(BACTOPIA_SAMPLESHEET.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
