/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules csvtk_concat as csvtk_concat_tsv, bracken as bracken_module, csvtk_concat as csvtk_concat_adjusted
 *
 * @input reads
 * Channel containing reads data
 *
 * @input database
 * Channel containing database data
 *
 * @output tsv                        Tsv
 * @output special_tsv                Special Tsv
 * @output classified                 Classified
 * @output unclassified               Unclassified
 * @output kraken2_report             Kraken2 Report
 * @output kraken2_output             Kraken2 Output
 * @output bracken_report             Bracken Report
 * @output krona                      Krona
 * @output abundances                 Abundances
 * @output classification             Classification
 * @output adjusted_abundances        Adjusted Abundances
 * @output merged_tsv                 Merged Tsv
 * @output merged_adjusted_abundances Merged Adjusted Abundances
 * @output results                    Aggregated results channel containing all output files
 * @output logs                       Aggregated logs channel containing all execution logs
 * @output nf_logs                    Aggregated Nextflow execution logs from all processes
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
    reads: Channel<Tuple<Map, Set<Path>>>
    database: Path

    main:
    BRACKEN_MODULE(reads, database)

    // Merge Bracken Primary/Secondary Species abundance
    CSVTK_CONCAT_TSV(gather(BRACKEN_MODULE.out.tsv, 'bracken-species-abundance'), 'tsv', 'tsv')

    // Merge Bracken adjusted abundance
    CSVTK_CONCAT_ADJUSTED(gather(BRACKEN_MODULE.out.adjusted_abundances, 'bracken-adjusted'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.tsv
    special_tsv: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.special_tsv
    classified: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.classified
    unclassified: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.unclassified
    kraken2_report: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.kraken2_report
    kraken2_output: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.kraken2_output
    bracken_report: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.bracken_report
    krona: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.krona
    abundances: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.abundances
    classification: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.classification
    adjusted_abundances: Channel<Tuple<Map, Path>> = BRACKEN_MODULE.out.adjusted_abundances
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_TSV.out.csv
    merged_adjusted_abundances: Channel<Tuple<Map, Path>> = CSVTK_CONCAT_ADJUSTED.out.csv

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
