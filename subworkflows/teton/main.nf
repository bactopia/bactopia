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
    // Remove host reads
    SCRUBBER(reads)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)
    ch_logs = ch_logs.mix(SCRUBBER.out.logs)
    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed)
    ch_versions = ch_versions.mix(BRACKEN.out.versions)
    ch_logs = ch_logs.mix(BRACKEN.out.logs)
    // Determine genome size and create sample sheet
    BACTOPIA_SAMPLESHEET(BRACKEN.out.classification)
    ch_versions = ch_versions.mix(BACTOPIA_SAMPLESHEET.out.versions)
    ch_logs = ch_logs.mix(BACTOPIA_SAMPLESHEET.out.logs)

    // Merge bacteria TSVs
    BACTOPIA_SAMPLESHEET.out.bacteria_tsv
        .collect{_meta, tsv -> tsv}
        .map{ tsv -> [[id:'teton-bacteria'], tsv]}
        .set{ ch_merge_bacteria }
    CSVTK_CONCAT(ch_merge_bacteria, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

    // Merge nonbacteria TSVs 
    BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv
        .collect{_meta, tsv -> tsv}
        .map{ tsv -> [[id:'teton-nonbacteria'], tsv]}
        .set{ ch_merge_nonbacteria }
    CSVTK_JOIN(ch_merge_nonbacteria, 'tsv', 'tsv', 'inner')
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    emit:
    bacteria_tsv = BACTOPIA_SAMPLESHEET.out.bacteria_tsv
    nonbacteria_tsv = BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv
    sizemeup = BACTOPIA_SAMPLESHEET.out.sizemeup
    logs = ch_logs
    nf_logs = BACTOPIA_SAMPLESHEET.out.nf_nf_begin.mix(
        BACTOPIA_SAMPLESHEET.out.nf_nf_err,
        BACTOPIA_SAMPLESHEET.out.nf_nf_log,
        BACTOPIA_SAMPLESHEET.out.nf_nf_out,
        BACTOPIA_SAMPLESHEET.out.nf_nf_run,
        BACTOPIA_SAMPLESHEET.out.nf_nf_sh,
        BACTOPIA_SAMPLESHEET.out.nf_nf_trace,
        BRACKEN.out.nf_logs,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        CSVTK_JOIN.out.nf_begin,
        CSVTK_JOIN.out.nf_err,
        CSVTK_JOIN.out.nf_log,
        CSVTK_JOIN.out.nf_out,
        CSVTK_JOIN.out.nf_run,
        CSVTK_JOIN.out.nf_sh,
        CSVTK_JOIN.out.nf_trace,
        SCRUBBER.out.nf_logs
    )
    versions = ch_versions
}
