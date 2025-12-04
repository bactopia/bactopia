//
// teton - Taxonomic classification and estimated species abundances
//
nextflow.preview.types = true

include { SCRUBBER                                 } from '../scrubber/main'
include { BRACKEN                                  } from '../bracken/main'
include { BACTOPIA_SAMPLESHEET                     } from '../../modules/bactopia/teton/main'
include { CSVTK_JOIN                               } from '../../modules/csvtk/join/main'
include { CSVTK_CONCAT                             } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_BACTERIA    } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_NONBACTERIA } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SIZEMEUP    } from '../../modules/csvtk/concat/main'
include { flattenPaths                             } from 'plugin/nf-bactopia'
include { gather                                   } from 'plugin/nf-bactopia'

workflow TETON {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>
    db: Path?
    use_srascrubber: Boolean

    main:
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
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
    ch_join_teton = SCRUBBER.out.special_tsv.join(BRACKEN.out.special_tsv, by:[0]).map{ meta, csv1, csv2 -> [meta, csv1, csv2] }
    CSVTK_JOIN(ch_join_teton, 'tsv', 'tsv', 'sample')
    ch_logs = ch_logs.mix(CSVTK_JOIN.out.logs)
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    // Merge reports
    CSVTK_CONCAT(gather(CSVTK_JOIN.out.csv, 'teton'), 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT.out.logs)
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    // Merge Teton prepare (bacteria)
    CSVTK_CONCAT_BACTERIA(gather(BACTOPIA_SAMPLESHEET.out.bacteria_tsv, 'teton-prepare'), 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT_BACTERIA.out.logs)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_BACTERIA.out.versions)

    // Merge Teton prepare (non-bacteria)
    CSVTK_CONCAT_NONBACTERIA(gather(BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv, 'teton-prepare-nonbacteria'), 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT_NONBACTERIA.out.logs)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_NONBACTERIA.out.versions)

    // Merge sizemeup results
    CSVTK_CONCAT_SIZEMEUP(gather(BACTOPIA_SAMPLESHEET.out.sizemeup, 'sizemeup'), 'tsv', 'tsv')
    ch_logs = ch_logs.mix(CSVTK_CONCAT_SIZEMEUP.out.logs)
    ch_versions = ch_versions.mix(CSVTK_CONCAT_SIZEMEUP.out.versions)

    emit:
    // Individual outputs
    bacteria_tsv: Channel<Tuple<Map, Path>> = BACTOPIA_SAMPLESHEET.out.bacteria_tsv
    merged_bacteria_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_BACTERIA.out.csv
    nonbacteria_tsv: Channel<Tuple<Map, Path>> = BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv
    merged_nonbacteria_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_NONBACTERIA.out.csv
    sizemeup: Channel<Tuple<Map, Path>> = BACTOPIA_SAMPLESHEET.out.sizemeup
    merged_sizemeup: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_SIZEMEUP.out.csv
    report: Channel<Tuple<Map, Path>> = CSVTK_JOIN.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_results,
        BACTOPIA_SAMPLESHEET.out.bacteria_tsv,
        BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv,
        BACTOPIA_SAMPLESHEET.out.sizemeup,
        CSVTK_CONCAT.out.csv,
        CSVTK_JOIN.out.csv,
        CSVTK_CONCAT_BACTERIA.out.csv,
        CSVTK_CONCAT_NONBACTERIA.out.csv,
        CSVTK_CONCAT_SIZEMEUP.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_nf_logs,
        BACTOPIA_SAMPLESHEET.out.nf_logs,
        CSVTK_JOIN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs,
        CSVTK_CONCAT_BACTERIA.out.nf_logs,
        CSVTK_CONCAT_NONBACTERIA.out.nf_logs,
        CSVTK_CONCAT_SIZEMEUP.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([ch_versions])
}
