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
 * @output tsv         LisSero serotype prediction results in TSV format
 * @output merged_tsv  Combined TSV file containing serotype results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { LISSERO as LISSERO_MODULE } from '../../modules/lissero/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow LISSERO {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    LISSERO_MODULE(assembly)
    CSVTK_CONCAT(gather(LISSERO_MODULE.out.tsv, 'lissero'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = LISSERO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        LISSERO_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        LISSERO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        LISSERO_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        LISSERO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
