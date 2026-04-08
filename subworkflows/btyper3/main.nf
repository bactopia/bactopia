/**
 * In silico taxonomic classification of Bacillus cereus group genomes.
 *
 * This subworkflow performs taxonomic classification of Bacillus cereus group
 * genomes using [BTyper3](https://github.com/lmc297/BTyper3), which provides
 * comprehensive classification including species, lineage, and toxin gene detection.
 * The results from individual samples are aggregated into a combined summary file.
 *
 * @status stable
 * @keywords Bacillus, cereus, taxonomy, typing, toxin genes
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation btyper3
 *
 * @modules csvtk_concat, btyper3
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for Bacillus cereus group classification
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited BTyper3 typing and classification results
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { BTYPER3 as BTYPER3_MODULE } from '../../modules/btyper3/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gatherCsvtk               } from 'plugin/nf-bactopia'

workflow BTYPER3 {
    take:
    assembly: Channel<Record>

    main:
    ch_btyper3 = BTYPER3_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_btyper3, 'tsv', [name: 'btyper3']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_btyper3
    run_outputs = ch_csvtk_concat
}
