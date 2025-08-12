//
// teton - Taxonomic classification and estimated species abundances
//
include { SCRUBBER } from '../scrubber/main'
include { BRACKEN } from '../bracken/main'
include { BACTOPIA_SAMPLESHEET } from '../../modules/bactopia/teton/main'
include { CSVTK_JOIN } from '../../modules/csvtk/join/main'
include { CSVTK_CONCAT } from '../../modules/csvtk/concat/main'

workflow TETON {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // Remove host reads
    SCRUBBER(reads)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)
    ch_logs = ch_logs.mix(SCRUBBER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SCRUBBER.out.nf_logs)

    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed)
    ch_versions = ch_versions.mix(BRACKEN.out.versions)
    ch_logs = ch_logs.mix(BRACKEN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BRACKEN.out.nf_logs)

    // Determine genome size and create sample sheet
    BACTOPIA_SAMPLESHEET(BRACKEN.out.classification)
    ch_versions = ch_versions.mix(BACTOPIA_SAMPLESHEET.out.versions)
    ch_logs = ch_logs.mix(BACTOPIA_SAMPLESHEET.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BACTOPIA_SAMPLESHEET.out.nf_logs)

    emit:
    bacteria_tsv = BACTOPIA_SAMPLESHEET.out.bacteria_tsv
    nonbacteria_tsv = BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv
    sizemeup = BACTOPIA_SAMPLESHEET.out.sizemeup
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
