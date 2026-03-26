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
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, database-dependent, conditional-logic
 * @citation kraken2, bracken
 *
 * @subworkflows scrubber, bracken
 * @modules bactopia_samplesheet, csvtk_join, csvtk_concat
 *
 * @input record(meta, r1, r2, se, lr)
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
 * @output scrubber_outputs
 * Per-sample scrubbing records from SCRUBBER
 *
 * @output bracken_outputs
 * Per-sample classification records from BRACKEN
 *
 * @output samplesheet_outputs
 * Per-sample genome size estimation and routing records
 *
 * @output report_outputs
 * Joined scrubber + classification report record from CSVTK_JOIN
 *
 * @output merged_report
 * Merged reports across all samples from CSVTK_CONCAT
 *
 * @output merged_bacteria
 * Merged bacterial sample sheets from CSVTK_CONCAT
 *
 * @output merged_nonbacteria
 * Merged non-bacterial sample sheets from CSVTK_CONCAT
 *
 * @output merged_sizemeup
 * Merged genome size estimates from CSVTK_CONCAT
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
include { gather                                   } from 'plugin/nf-bactopia'

workflow TETON {
    take:
    reads: Channel<Record>
    db: Path?
    use_srascrubber: Boolean
    nohuman_db: Path?
    download_nohuman: Boolean

    main:
    // Remove host reads
    SCRUBBER(reads, use_srascrubber, nohuman_db, download_nohuman)

    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed, db)

    // Determine genome size and create sample sheet
    ch_classification = BRACKEN.out.sample_outputs.map { r -> record(_meta: r.meta, classification: r.classification) }
    BACTOPIA_SAMPLESHEET(ch_classification)

    // Join Scrubber and Bracken results
    ch_bracken_special = BRACKEN.out.sample_outputs.map { r -> record(special_meta: r.special_meta, tsv: r.tsv) }
    ch_join_teton = SCRUBBER.out.special_tsv.join(ch_bracken_special, by:[0]).map { meta, csv1, csv2 -> record(meta: meta, csv1: csv1, csv2: csv2) }
    CSVTK_JOIN(ch_join_teton, 'tsv', 'tsv', 'sample')

    // Merge reports
    CSVTK_CONCAT(gather(CSVTK_JOIN.out, 'csv', [name: 'teton']), 'tsv', 'tsv')

    // Merge Teton prepare (bacteria)
    CSVTK_CONCAT_BACTERIA(gather(BACTOPIA_SAMPLESHEET.out, 'bacteria_tsv', [name: 'teton-prepare']), 'tsv', 'tsv')

    // Merge Teton prepare (non-bacteria)
    CSVTK_CONCAT_NONBACTERIA(gather(BACTOPIA_SAMPLESHEET.out, 'nonbacteria_tsv', [name: 'teton-prepare-nonbacteria']), 'tsv', 'tsv')

    // Merge sizemeup results
    CSVTK_CONCAT_SIZEMEUP(gather(BACTOPIA_SAMPLESHEET.out, 'sizemeup', [name: 'sizemeup']), 'tsv', 'tsv')

    emit:
    // Published outputs
    scrubber_outputs = SCRUBBER.out.sample_outputs
    bracken_outputs = BRACKEN.out.sample_outputs
    samplesheet_outputs = BACTOPIA_SAMPLESHEET.out
    report_outputs = CSVTK_JOIN.out
    merged_report = CSVTK_CONCAT.out
    merged_bacteria = CSVTK_CONCAT_BACTERIA.out
    merged_nonbacteria = CSVTK_CONCAT_NONBACTERIA.out
    merged_sizemeup = CSVTK_CONCAT_SIZEMEUP.out
}
