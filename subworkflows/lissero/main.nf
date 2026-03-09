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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for L. monocytogenes serotype prediction
 *
 * @output sample_outputs  Per-sample records containing meta, tsv, results, logs, nf_logs, versions
 * @output run_outputs   Cross-sample aggregation record
 */
nextflow.preview.types = true

include { LISSERO as LISSERO_MODULE } from '../../modules/lissero/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { gather                    } from 'plugin/nf-bactopia'

workflow LISSERO {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    LISSERO_MODULE(assembly)
    CSVTK_CONCAT(gather(LISSERO_MODULE.out, 'lissero', field: 'tsv'), 'tsv', 'tsv')

    emit:
    // Per-sample records (contains meta, tsv, results, logs, nf_logs, versions)
    sample_outputs = LISSERO_MODULE.out
    // Cross-sample aggregation record
    run_outputs = CSVTK_CONCAT.out
}
