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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for H. parasuis serotype prediction
 *
 * @output sample_outputs  Per-sample records from HPSUISSERO_MODULE
 * @output run_outputs   Cross-sample aggregation record from CSVTK_CONCAT
 */
nextflow.preview.types = true

include { HPSUISSERO as HPSUISSERO_MODULE } from '../../modules/hpsuissero/main'
include { CSVTK_CONCAT                    } from '../../modules/csvtk/concat/main'
include { gather                          } from 'plugin/nf-bactopia'

workflow HPSUISSERO {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    HPSUISSERO_MODULE(assembly)
    CSVTK_CONCAT(gather(HPSUISSERO_MODULE.out, 'hpsuissero', field: 'tsv'), 'tsv', 'tsv')

    emit:
    // Per-sample records
    sample_outputs = HPSUISSERO_MODULE.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
