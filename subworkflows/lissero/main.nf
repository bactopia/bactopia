/**
 * In silico serotype prediction for Listeria monocytogenes.
 *
 * This subworkflow performs serotype prediction for Listeria monocytogenes
 * using [LisSero](https://github.com/MDU-PHL/LisSero), which identifies specific
 * serotype markers in genome assemblies. The tool provides rapid classification
 * into the major L. monocytogenes serotypes, which is important for outbreak
 * investigation and tracking.
 *
 * @status stable
 * @keywords Listeria, monocytogenes, serotype, outbreak
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation lissero
 *
 * @modules csvtk_concat, lissero
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for L. monocytogenes serotype prediction
 *
 * @output sample_outputs
 * - `tsv`: Tab-delimited LisSero results with predicted serogroup and marker gene detection
 *
 * @output run_outputs
 * - `csv`: Aggregated results in CSV format
 */
nextflow.preview.types = true

include { LISSERO as LISSERO_MODULE } from '../../modules/lissero/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gatherCsvtk                    } from 'plugin/nf-bactopia'

workflow LISSERO {
    take:
    assembly: Channel<Record>

    main:
    LISSERO_MODULE(assembly)
    CSVTK_CONCAT(gatherCsvtk(LISSERO_MODULE.out, 'tsv', [name: 'lissero']), 'tsv', 'tsv')

    emit:
    // Published outputs
    sample_outputs = LISSERO_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
