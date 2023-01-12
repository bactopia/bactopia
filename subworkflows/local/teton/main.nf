//
// teton - Taxonomic classification and estimated species abundances
//
include { SCRUBBER } from '../scrubber/main'
include { KRAKEN2_BRACKEN } from '../bracken/main'
include { MIDAS } from '../midas/main'
include { CSVTK_JOIN } from '../../../modules/nf-core/csvtk/join/main' addParams( options: [publish_to_base: true, logs_subdir: 'amrfinderplus-proteins'] )
include { CSVTK_CONCAT } from '../../../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true, logs_subdir: 'amrfinderplus-proteins'] )

workflow TETON {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    if (params.skip_scrubber) {
        // Taxon Classification & Abundance
        KRAKEN2_BRACKEN(reads)
        ch_versions = ch_versions.mix(KRAKEN2_BRACKEN.out.versions)

        if (!params.skip_midas) {
            // Species Abundance
            MIDAS(reads)
            ch_versions = ch_versions.mix(MIDAS.out.versions)
        }
    } else {
        // Remove host reads
        SCRUBBER(reads)
        ch_versions = ch_versions.mix(SCRUBBER.out.versions)

        // Taxon Classification & Abundance
        KRAKEN2_BRACKEN(SCRUBBER.out.scrubbed)
        ch_versions = ch_versions.mix(KRAKEN2_BRACKEN.out.versions)

        if (!params.skip_midas) {
            // Species Abundance
            MIDAS(SCRUBBER.out.scrubbed)
            ch_versions = ch_versions.mix(MIDAS.out.versions)
        }
    }

    // Summarize the results

    emit:

    versions = ch_versions // channel: [ versions.yml ]
}
