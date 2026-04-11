/**
 * Perform taxonomic classification and estimate bacterial genome sizes.
 *
 * This subworkflow processes raw sequencing reads through a taxonomic classification
 * pipeline using [Kraken2](https://github.com/DerrickWood/kraken2) and [Bracken](https://github.com/jenniferlu717/Bracken)
 * to estimate bacterial genome sizes and separate bacterial from non-bacterial organisms.
 * It first removes host reads using the scrubber subworkflow, then classifies reads,
 * and finally creates sample sheets with genome size estimates for downstream Bactopia analysis.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomy, classification, kraken, bracken, genome size
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation,database-dependent,conditional-logic
 * @citation kraken2, bracken
 *
 * @modules bactopia_teton, csvtk_join, csvtk_concat
 * @subworkflows scrubber, bracken
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input db
 * Optional Kraken2 database path for taxonomic classification
 *
 * @input use_srascrubber
 * Boolean flag to use SRA scrubber for host read removal
 *
 * @output sample_outputs
 *
 * @output run_outputs
 */
// bactopia-lint: ignore S013,S015
nextflow.preview.types = true

include { SCRUBBER                                 } from '../scrubber/main'
include { BRACKEN                                  } from '../bracken/main'
include { BACTOPIA_SAMPLESHEET                     } from '../../modules/bactopia/teton/main'
include { CSVTK_JOIN                               } from '../../modules/csvtk/join/main'
include { CSVTK_CONCAT                             } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_BACTERIA    } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_NONBACTERIA } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_SIZEMEUP    } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                              } from 'plugin/nf-bactopia'

workflow TETON {
    take:
    reads: Channel<Record>
    db: Value<Path?>
    use_srascrubber: Boolean
    nohuman_db: Value<Path?>
    download_nohuman: Boolean
    nohuman_save_as_tarball: Boolean

    main:
    // Remove host reads
    ch_scrubber = SCRUBBER(reads, use_srascrubber, nohuman_db, download_nohuman, nohuman_save_as_tarball)

    // Taxon Classification & Abundance
    ch_bracken = BRACKEN(ch_scrubber.scrubbed, db)

    // Determine genome size and create sample sheet
    ch_classification = ch_bracken.sample_outputs.map { r -> record(meta: r.meta, classification: r.classification) }
    ch_bactopia_samplesheet = BACTOPIA_SAMPLESHEET(ch_classification)

    // Join Scrubber and Bracken results
    ch_bracken_special = ch_bracken.sample_outputs.map { r -> record(special_meta: r.special_meta, tsv: r.tsv) }
    ch_join_teton = ch_scrubber.special_tsv.join(ch_bracken_special, by: 'special_meta').map { r -> record(meta: r.special_meta, csv1: r.scrub_report, csv2: r.tsv) }
    ch_csvtk_join = CSVTK_JOIN(ch_join_teton, 'tsv', 'tsv', 'sample')

    // Merge reports
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_csvtk_join, 'csv', [name: 'teton']), 'tsv', 'tsv')

    // Merge Teton prepare (bacteria)
    ch_csvtk_concat_bacteria = CSVTK_CONCAT_BACTERIA(gatherCsvtk(ch_bactopia_samplesheet, 'bacteria_tsv', [name: 'teton-prepare']), 'tsv', 'tsv')

    // Merge Teton prepare (non-bacteria)
    ch_csvtk_concat_nonbacteria = CSVTK_CONCAT_NONBACTERIA(gatherCsvtk(ch_bactopia_samplesheet, 'nonbacteria_tsv', [name: 'teton-prepare-nonbacteria']), 'tsv', 'tsv')

    // Merge sizemeup results
    ch_csvtk_concat_sizemeup = CSVTK_CONCAT_SIZEMEUP(gatherCsvtk(ch_bactopia_samplesheet, 'sizemeup', [name: 'sizemeup']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_scrubber.sample_outputs.mix(
        ch_bracken.sample_outputs,
        ch_bactopia_samplesheet,
        ch_csvtk_join
    )
    run_outputs = ch_scrubber.run_outputs.mix(
        ch_bracken.run_outputs,
        ch_csvtk_concat,
        ch_csvtk_concat_bacteria,
        ch_csvtk_concat_nonbacteria,
        ch_csvtk_concat_sizemeup
    )
}
