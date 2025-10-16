//
// teton - Taxonomic classification and estimated species abundances
//
include { SCRUBBER             } from '../scrubber/main'
include { BRACKEN              } from '../bracken/main'
include { BACTOPIA_SAMPLESHEET } from '../../modules/bactopia/teton/main'
include { CSVTK_JOIN           } from '../../modules/csvtk/join/main'
include { CSVTK_CONCAT         } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_BACTERIA    } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_NONBACTERIA } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SIZEMEUP    } from '../../modules/csvtk/concat/main'

workflow TETON {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    db
    use_srascrubber

    main:
    ch_results = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_versions = Channel.empty()

    // Remove host reads
    SCRUBBER(reads, use_srascrubber)
    ch_results = ch_results.mix(SCRUBBER.out.results)
    ch_logs = ch_logs.mix(SCRUBBER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SCRUBBER.out.nf_logs)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)

    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed, db)
    ch_results = ch_results.mix(BRACKEN.out.results)
    ch_logs = ch_logs.mix(BRACKEN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(BRACKEN.out.nf_logs)
    ch_versions = ch_versions.mix(BRACKEN.out.versions)

    // Determine genome size and create sample sheet
    BACTOPIA_SAMPLESHEET(BRACKEN.out.classification)
    ch_logs = ch_logs.mix(BACTOPIA_SAMPLESHEET.out.logs)
    ch_versions = ch_versions.mix(BACTOPIA_SAMPLESHEET.out.versions)

    // Join Scrubber and Bracken results
    SCRUBBER.out.special_tsv.join(BRACKEN.out.special_tsv, by:[0]).map{ meta, csv1, csv2 -> [meta, csv1, csv2] }.set{ ch_join_teton }
    CSVTK_JOIN(ch_join_teton, 'tsv', 'tsv', 'sample')
    ch_logs = ch_logs.mix(CSVTK_JOIN.out.logs)
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    // Merge reports
    CSVTK_JOIN.out.csv.collect{_meta, csv -> csv}.map{ csv -> [[id:'teton'], csv]}.set{ ch_merge_teton }
    CSVTK_CONCAT(ch_merge_teton, 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    // Merge Teton prepare (bacteria)
    BACTOPIA_SAMPLESHEET.out.bacteria_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'teton-prepare'], tsv]}.set{ ch_merge_prepare }
    CSVTK_CONCAT_BACTERIA(ch_merge_prepare, 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT_BACTERIA.out.logs)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_BACTERIA.out.versions)

    // Merge Teton prepare (non-bacteria)
    BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'teton-prepare-nonbacteria'], tsv]}.set{ ch_merge_prepare_non }
    CSVTK_CONCAT_NONBACTERIA(ch_merge_prepare_non, 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT_NONBACTERIA.out.logs)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_NONBACTERIA.out.versions)

    // Merge sizemeup results
    BACTOPIA_SAMPLESHEET.out.sizemeup.collect{_meta, tsv -> tsv}.map{ tsv -> [[id:'sizemeup'], tsv]}.set{ ch_merge_sizemeup }
    CSVTK_CONCAT_SIZEMEUP(ch_merge_sizemeup, 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT_SIZEMEUP.out.logs)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SIZEMEUP.out.versions)

    emit:
    // Individual outputs
    bacteria_tsv = BACTOPIA_SAMPLESHEET.out.bacteria_tsv
    merged_bacteria_tsv = CSVTK_CONCAT_BACTERIA.out.csv
    nonbacteria_tsv = BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv
    merged_nonbacteria_tsv = CSVTK_CONCAT_NONBACTERIA.out.csv
    sizemeup = BACTOPIA_SAMPLESHEET.out.sizemeup
    merged_sizemeup = CSVTK_CONCAT_SIZEMEUP.out.csv
    report = CSVTK_JOIN.out.csv

    // Generic aggregate outputs
    results = ch_results.mix(
        BACTOPIA_SAMPLESHEET.out.bacteria_tsv,
        BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv,
        BACTOPIA_SAMPLESHEET.out.sizemeup,
        CSVTK_CONCAT.out.csv,
        CSVTK_JOIN.out.csv,
        CSVTK_CONCAT_BACTERIA.out.csv,
        CSVTK_CONCAT_NONBACTERIA.out.csv,
        CSVTK_CONCAT_SIZEMEUP.out.csv
    )
    logs = ch_logs
    nf_logs = ch_nf_logs.mix(
        BACTOPIA_SAMPLESHEET.out.nf_begin,
        BACTOPIA_SAMPLESHEET.out.nf_err,
        BACTOPIA_SAMPLESHEET.out.nf_log,
        BACTOPIA_SAMPLESHEET.out.nf_out,
        BACTOPIA_SAMPLESHEET.out.nf_run,
        BACTOPIA_SAMPLESHEET.out.nf_sh,
        BACTOPIA_SAMPLESHEET.out.nf_trace,
        CSVTK_JOIN.out.nf_begin,
        CSVTK_JOIN.out.nf_err,
        CSVTK_JOIN.out.nf_log,
        CSVTK_JOIN.out.nf_out,
        CSVTK_JOIN.out.nf_run,
        CSVTK_JOIN.out.nf_sh,
        CSVTK_JOIN.out.nf_trace,
        CSVTK_CONCAT.out.nf_begin,
        CSVTK_CONCAT.out.nf_err,
        CSVTK_CONCAT.out.nf_log,
        CSVTK_CONCAT.out.nf_out,
        CSVTK_CONCAT.out.nf_run,
        CSVTK_CONCAT.out.nf_sh,
        CSVTK_CONCAT.out.nf_trace,
        CSVTK_CONCAT_BACTERIA.out.nf_begin,
        CSVTK_CONCAT_BACTERIA.out.nf_err,
        CSVTK_CONCAT_BACTERIA.out.nf_log,
        CSVTK_CONCAT_BACTERIA.out.nf_out,
        CSVTK_CONCAT_BACTERIA.out.nf_run,
        CSVTK_CONCAT_BACTERIA.out.nf_sh,
        CSVTK_CONCAT_BACTERIA.out.nf_trace,
        CSVTK_CONCAT_NONBACTERIA.out.nf_begin,
        CSVTK_CONCAT_NONBACTERIA.out.nf_err,
        CSVTK_CONCAT_NONBACTERIA.out.nf_log,
        CSVTK_CONCAT_NONBACTERIA.out.nf_out,
        CSVTK_CONCAT_NONBACTERIA.out.nf_run,
        CSVTK_CONCAT_NONBACTERIA.out.nf_sh,
        CSVTK_CONCAT_NONBACTERIA.out.nf_trace,
        CSVTK_CONCAT_SIZEMEUP.out.nf_begin,
        CSVTK_CONCAT_SIZEMEUP.out.nf_err,
        CSVTK_CONCAT_SIZEMEUP.out.nf_log,
        CSVTK_CONCAT_SIZEMEUP.out.nf_out,
        CSVTK_CONCAT_SIZEMEUP.out.nf_run,
        CSVTK_CONCAT_SIZEMEUP.out.nf_sh,
        CSVTK_CONCAT_SIZEMEUP.out.nf_trace
    )
    versions = ch_versions
}
