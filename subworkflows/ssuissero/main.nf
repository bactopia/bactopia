/**
 * Predict serotypes of Streptococcus suis from genome assemblies.
 *
 * This subworkflow uses [SsuisSero](https://github.com/jimmyliu1326/SsuisSero) to predict
 * serotypes of *Streptococcus suis* strains from genome assemblies based on the presence
 * of specific capsular genes. It processes each sample individually and aggregates the
 * results into a single consolidated report.
 *
 * @status stable
 * @keywords streptococcus suis, serotype, typing, prediction, capsular genes
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation ssuissero
 *
 * @modules ssuissero, csvtk_concat
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output sample_outputs
 * - `tsv`: SsuisSero results in TSV format
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { SSUISSERO as SSUISSERO_MODULE } from '../../modules/ssuissero/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                   } from 'plugin/nf-bactopia'

workflow SSUISSERO {
    take:
    assembly: Channel<Record>

    main:
    ch_ssuissero = SSUISSERO_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_ssuissero, 'tsv', [name: 'ssuissero']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_ssuissero
    run_outputs = ch_csvtk_concat
}
