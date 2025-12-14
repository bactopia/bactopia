/**
 * In silico Legionella pneumophila Sequence Based Typing.
 *
 * This subworkflow performs sequence-based typing of Legionella pneumophila
 * using [legsta](https://github.com/tseemann/legsta), which identifies the
 * Sequence Type (ST) based on the seven-locus scheme. The tool analyzes
 * allele profiles and provides epidemiological typing data for outbreak
 * investigation and population studies.
 *
 * @status stable
 * @keywords Legionella, pneumophila, sequence typing, ST, epidemiology
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation legsta
 *
 * @modules csvtk_concat, legsta
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for L. pneumophila sequence typing
 *
 * @output tsv         legsta sequence typing results with ST assignments
 * @output merged_tsv  Combined TSV file containing typing results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { LEGSTA as LEGSTA_MODULE } from '../../modules/legsta/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow LEGSTA {
    take:
    assembly: Channel<Tuple<Map, Path>>

    main:
    LEGSTA_MODULE(assembly)
    CSVTK_CONCAT(gather(LEGSTA_MODULE.out.tsv, 'legsta'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = LEGSTA_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        LEGSTA_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        LEGSTA_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        LEGSTA_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        LEGSTA_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
