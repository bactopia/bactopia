/**
 * Estimate species abundance from metagenomic reads.
 *
 * This subworkflow performs taxonomic classification and abundance estimation using [Kraken2](https://github.com/DerrickWood/kraken2)
 * and [Bracken](https://github.com/jenniferlu717/Bracken). It processes metagenomic reads, classifies them against a reference database,
 * and generates abundance estimates at different taxonomic levels with optional abundance correction.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomic classification, abundance estimation, kraken2, bracken
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, conditional-logic, aggregation
 * @citation kraken2, bracken
 *
 * @modules bracken, csvtk_concat as csvtk_concat_tsv, csvtk_concat as csvtk_concat_adjusted
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input database
 * Path to the Kraken2 database for taxonomic classification.
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited summary of Bracken primary and secondary species abundances
 * - `special_meta`: A simplified metadata map for internal use
 * - `classified`: Reads classified to belong to any of the taxa on the Kraken2 database
 * - `unclassified`: Reads not classified to belong to any of the taxa on the Kraken2 database
 * - `kraken2_report`: Kraken2 report containing stats about classified and not classified reads
 * - `kraken2_output`: Kraken2 output file containing the taxonomic classification of each read
 * - `bracken_report`: Bracken report containing stats about classified and not classified reads
 * - `krona`: Interactive Krona HTML visualization
 * - `abundances`: Bracken abundance estimates for each taxon
 * - `classification`: Bracken per-read classification details
 * - `adjusted_abundances`: Bracken abundance estimates adjusted for unclassified reads
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BRACKEN as BRACKEN_MODULE             } from '../../modules/bracken/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_TSV      } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_ADJUSTED } from '../../modules/csvtk/concat/main'
include { gather                                } from 'plugin/nf-bactopia'

workflow BRACKEN {
    take:
    reads: Channel<Record>
    database: Path

    main:
    BRACKEN_MODULE(reads, database)

    // Merge Bracken Primary/Secondary Species abundance
    CSVTK_CONCAT_TSV(gather(BRACKEN_MODULE.out, 'tsv', [name: 'bracken-species-abundance']), 'tsv', 'tsv')

    // Merge Bracken adjusted abundance
    CSVTK_CONCAT_ADJUSTED(gather(BRACKEN_MODULE.out, 'adjusted_abundances', [name: 'bracken-adjusted']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = BRACKEN_MODULE.out
    run_outputs = CSVTK_CONCAT_TSV.out.mix(CSVTK_CONCAT_ADJUSTED.out)
}
