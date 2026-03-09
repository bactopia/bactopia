/**
 * Perform taxonomic classification and estimate bacterial genome sizes.
 *
 * This subworkflow processes raw sequencing reads through a taxonomic classification
 * pipeline using [Kraken2](https://github.com/DerrickWood/kraken2) and [Bracken](https://github.com/jenniferlu717/Bracken)
 * to estimate bacterial genome sizes and separate bacterial from non-bacterial organisms.
 * It first removes host reads using the scrubber subworkflow, then classifies reads,
 * and finally creates sample sheets with genome size estimates for downstream Bactopia analysis.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomy, classification, kraken, bracken, genome size
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, database-dependent, conditional-logic
 * @citation kraken2, bracken
 *
 * @subworkflows scrubber, bracken
 * @modules bactopia_samplesheet, csvtk_join, csvtk_concat
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input db
 * Optional Kraken2 database path for taxonomic classification
 *
 * @input use_srascrubber
 * Boolean flag to use SRA scrubber for host read removal
 *
 * @output bacteria_tsv           Per-sample TSV files containing bacterial organisms and their properties
 * @output merged_bacteria_tsv    Consolidated TSV file of all bacterial organisms across samples
 * @output nonbacteria_tsv        Per-sample TSV files containing non-bacterial organisms
 * @output merged_nonbacteria_tsv Consolidated TSV file of all non-bacterial organisms across samples
 * @output sizemeup               Per-sample TSV files with genome size estimates
 * @output merged_sizemeup        Consolidated TSV file of genome size estimates across samples
 * @output report                 Joined TSV file combining scrubber and classification results
 * @output results                Aggregated results channel containing all output files
 * @output logs                   Aggregated logs channel containing all execution logs
 * @output nf_logs                Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions               Aggregated version information from all executed tools
 */
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
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    db: Path?
    use_srascrubber: Boolean

    main:
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

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
    ch_concat_logs = CSVTK_CONCAT.out.map { r -> tuple(r.meta, r.logs) }
    ch_logs = ch_logs.mix(ch_concat_logs)
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    // Merge Teton prepare (bacteria)
    CSVTK_CONCAT_BACTERIA(gather(BACTOPIA_SAMPLESHEET.out.bacteria_tsv, 'teton-prepare'), 'tsv', 'tsv')
    ch_bacteria_logs = CSVTK_CONCAT_BACTERIA.out.map { r -> tuple(r.meta, r.logs) }
    ch_bacteria_versions = CSVTK_CONCAT_BACTERIA.out.map { r -> tuple(r.meta, r.versions) }
    ch_logs = ch_logs.mix(ch_bacteria_logs)
    ch_versions = ch_versions.mix(ch_bacteria_versions)

    // Merge Teton prepare (non-bacteria)
    CSVTK_CONCAT_NONBACTERIA(gather(BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv, 'teton-prepare-nonbacteria'), 'tsv', 'tsv')
    ch_nonbacteria_logs = CSVTK_CONCAT_NONBACTERIA.out.map { r -> tuple(r.meta, r.logs) }
    ch_nonbacteria_versions = CSVTK_CONCAT_NONBACTERIA.out.map { r -> tuple(r.meta, r.versions) }
    ch_logs = ch_logs.mix(ch_nonbacteria_logs)
    ch_versions = ch_versions.mix(ch_nonbacteria_versions)

    // Merge sizemeup results
    CSVTK_CONCAT_SIZEMEUP(gather(BACTOPIA_SAMPLESHEET.out.sizemeup, 'sizemeup'), 'tsv', 'tsv')
    ch_sizemeup_logs = CSVTK_CONCAT_SIZEMEUP.out.map { r -> tuple(r.meta, r.logs) }
    ch_sizemeup_versions = CSVTK_CONCAT_SIZEMEUP.out.map { r -> tuple(r.meta, r.versions) }
    ch_logs = ch_logs.mix(ch_sizemeup_logs)
    ch_versions = ch_versions.mix(ch_sizemeup_versions)

    // Extract tuple channels from CSVTK_CONCAT record outputs for flattenPaths/emit compatibility
    ch_concat_csv = CSVTK_CONCAT.out.map { r -> tuple(r.meta, r.csv) }
    ch_concat_nf_logs = CSVTK_CONCAT.out.map { r -> tuple(r.meta, r.nf_logs) }
    ch_bacteria_csv = CSVTK_CONCAT_BACTERIA.out.map { r -> tuple(r.meta, r.csv) }
    ch_bacteria_nf_logs = CSVTK_CONCAT_BACTERIA.out.map { r -> tuple(r.meta, r.nf_logs) }
    ch_nonbacteria_csv = CSVTK_CONCAT_NONBACTERIA.out.map { r -> tuple(r.meta, r.csv) }
    ch_nonbacteria_nf_logs = CSVTK_CONCAT_NONBACTERIA.out.map { r -> tuple(r.meta, r.nf_logs) }
    ch_sizemeup_csv = CSVTK_CONCAT_SIZEMEUP.out.map { r -> tuple(r.meta, r.csv) }
    ch_sizemeup_nf_logs = CSVTK_CONCAT_SIZEMEUP.out.map { r -> tuple(r.meta, r.nf_logs) }

    emit:
    // Individual outputs
    bacteria_tsv: Channel<Tuple<Map, Set<Path>>> = BACTOPIA_SAMPLESHEET.out.bacteria_tsv
    merged_bacteria_tsv = ch_bacteria_csv
    nonbacteria_tsv: Channel<Tuple<Map, Set<Path>>> = BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv
    merged_nonbacteria_tsv = ch_nonbacteria_csv
    sizemeup: Channel<Tuple<Map, Set<Path>>> = BACTOPIA_SAMPLESHEET.out.sizemeup
    merged_sizemeup = ch_sizemeup_csv
    report: Channel<Tuple<Map, Set<Path>>> = CSVTK_JOIN.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_results,
        BACTOPIA_SAMPLESHEET.out.bacteria_tsv,
        BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv,
        BACTOPIA_SAMPLESHEET.out.sizemeup,
        ch_concat_csv,
        CSVTK_JOIN.out.csv,
        ch_bacteria_csv,
        ch_nonbacteria_csv,
        ch_sizemeup_csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_nf_logs,
        BACTOPIA_SAMPLESHEET.out.nf_logs,
        CSVTK_JOIN.out.nf_logs,
        ch_concat_nf_logs,
        ch_bacteria_nf_logs,
        ch_nonbacteria_nf_logs,
        ch_sizemeup_nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([ch_versions])
}
