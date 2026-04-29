/**
 * Scripts for finding and processing promoter variants upstream of mcr-1.
 *
 * This subworkflow identifies and characterizes promoter variants upstream of
 * the mcr-1 colistin resistance gene using [mcroni](https://github.com/liampshaw/mcroni).
 * The tool searches for mutations in the promoter region that may affect expression
 * levels of mcr-1, which is important for understanding the regulation of
 * plasmid-mediated colistin resistance.
 *
 * @status stable
 * @keywords mcr-1, colistin, resistance, promoter, variant
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation mcroni
 *
 * @modules csvtk_concat, mcroni
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembly files in FASTA format for mcr-1 promoter analysis
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited mcr-1 gene variation results
 * - `fa`: Extracted mcr-1 gene sequence in FASTA format (optional)
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.enable.types = true

include { MCRONI as MCRONI_MODULE } from '../../modules/mcroni/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { gatherCsvtk             } from 'plugin/nf-bactopia'

workflow MCRONI {
    take:
    assembly: Channel<Record>

    main:
    ch_mcroni = MCRONI_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_mcroni, 'tsv', [name: 'mcroni']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_mcroni
    run_outputs = ch_csvtk_concat
}
