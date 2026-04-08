/**
 * Rapid Haemophilus parasuis serotyping.
 *
 * This subworkflow performs serotyping of Haemophilus parasuis using
 * [HpsuisSero](https://github.com/jimmyliu1326/HpsuisSero), which identifies
 * serotype-specific markers in genome assemblies. The tool provides rapid
 * classification of H. parasuis isolates into their respective serotypes,
 * which is important for epidemiological surveillance and vaccine development.
 *
 * @status stable
 * @keywords Haemophilus, parasuis, serotype, epidemiology
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation hpsuissero
 *
 * @modules csvtk_concat, hpsuissero
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for H. parasuis serotype prediction
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited HpsuisSero serotype prediction results
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                     } from 'plugin/nf-bactopia'

workflow HPSUISSERO {
    take:
    assembly: Channel<Record>

    main:
    ch_hpsuissero = HPSUISSERO_MODULE(assembly)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_hpsuissero, 'tsv', [name: 'hpsuissero']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = ch_hpsuissero
    run_outputs = ch_csvtk_concat
}
