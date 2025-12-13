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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing serotype predictions
 * @output merged_tsv  Consolidated TSV file containing serotype predictions from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SSUISSERO as SSUISSERO_MODULE } from '../../modules/ssuissero/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow SSUISSERO {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    SSUISSERO_MODULE(assembly)
    CSVTK_CONCAT(gather(SSUISSERO_MODULE.out.tsv, 'ssuissero'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = SSUISSERO_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SSUISSERO_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
