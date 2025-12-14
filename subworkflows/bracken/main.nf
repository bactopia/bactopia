/**
 * Estimate species abundance from metagenomic reads.
 *
 * This subworkflow performs taxonomic classification and abundance estimation using [Kraken2](https://github.com/DerrickWood/kraken2)
 * and [Bracken](https://github.com/jenniferlu717/Bracken). It processes metagenomic reads, classifies them against a reference database,
 * and generates abundance estimates at different taxonomic levels with optional abundance correction.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomic classification, abundance estimation, kraken2, bracken
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, conditional-logic, aggregation
 * @citation kraken2, bracken
 *
 * @modules bracken, csvtk_concat as csvtk_concat_tsv, csvtk_concat as csvtk_concat_adjusted
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @input database
 * Path to the Kraken2 database for taxonomic classification.
 *
 * @output classified                 Taxonomically classified reads in FASTA format
 * @output unclassified               Unclassified reads not assigned to any taxa
 * @output kraken2_report             Kraken2 classification report with read counts per taxon
 * @output kraken2_output             Detailed Kraken2 classification output for each read
 * @output bracken_report             Bracken abundance estimation report
 * @output krona                      Krona-compatible XML file for interactive taxonomic visualization
 * @output abundances                 Species-level abundance estimates from Bracken
 * @output adjusted_abundances        Read count-adjusted abundance estimates
 * @output merged_tsv                 Merged abundance TSV files across all samples
 * @output merged_adjusted_abundances Merged adjusted abundance TSV files across all samples
 * @output results                    Aggregated results channel containing all output files
 * @output logs                       Aggregated logs channel containing all execution logs
 * @output nf_logs                    Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions                   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BRACKEN as BRACKEN_MODULE             } from '../../modules/bracken/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_TSV      } from '../../modules/csvtk/concat/main'
include { CSVTK_CONCAT as CSVTK_CONCAT_ADJUSTED } from '../../modules/csvtk/concat/main'
include { flattenPaths                          } from 'plugin/nf-bactopia'
include { gather                                } from 'plugin/nf-bactopia'

workflow BRACKEN {
    take:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    database: Path

    main:
    BRACKEN_MODULE(reads, database)

    // Merge Bracken Primary/Secondary Species abundance
    CSVTK_CONCAT_TSV(gather(BRACKEN_MODULE.out.tsv, 'bracken-species-abundance'), 'tsv', 'tsv')

    // Merge Bracken adjusted abundance
    CSVTK_CONCAT_ADJUSTED(gather(BRACKEN_MODULE.out.adjusted_abundances, 'bracken-adjusted'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.tsv
    special_tsv: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.special_tsv
    classified: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.classified
    unclassified: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.unclassified
    kraken2_report: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.kraken2_report
    kraken2_output: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.kraken2_output
    bracken_report: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.bracken_report
    krona: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.krona
    abundances: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.abundances
    classification: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.classification
    adjusted_abundances: Channel<Tuple<Map, Set<Path>>> = BRACKEN_MODULE.out.adjusted_abundances
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT_TSV.out.csv
    merged_adjusted_abundances: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT_ADJUSTED.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BRACKEN_MODULE.out.tsv,
        BRACKEN_MODULE.out.classified,
        BRACKEN_MODULE.out.unclassified,
        BRACKEN_MODULE.out.kraken2_report,
        BRACKEN_MODULE.out.kraken2_output,
        BRACKEN_MODULE.out.bracken_report,
        BRACKEN_MODULE.out.krona,
        BRACKEN_MODULE.out.abundances,
        BRACKEN_MODULE.out.classification,
        BRACKEN_MODULE.out.adjusted_abundances,
        CSVTK_CONCAT_TSV.out.csv,
        CSVTK_CONCAT_ADJUSTED.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BRACKEN_MODULE.out.logs,
        CSVTK_CONCAT_TSV.out.logs,
        CSVTK_CONCAT_ADJUSTED.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BRACKEN_MODULE.out.nf_logs,
        CSVTK_CONCAT_TSV.out.nf_logs,
        CSVTK_CONCAT_ADJUSTED.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BRACKEN_MODULE.out.versions,
        CSVTK_CONCAT_TSV.out.versions,
        CSVTK_CONCAT_ADJUSTED.out.versions
    ])
}
