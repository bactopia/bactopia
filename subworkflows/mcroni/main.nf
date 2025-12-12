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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for mcr-1 promoter analysis
 *
 * @output tsv         mcroni analysis results with identified promoter variants
 * @output merged_tsv  Combined TSV file containing promoter variant results from all samples
 * @output fa          Extracted promoter sequences in FASTA format
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { MCRONI as MCRONI_MODULE } from '../../modules/mcroni/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow MCRONI {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    MCRONI_MODULE(assembly)
    CSVTK_CONCAT(gather(MCRONI_MODULE.out.tsv, 'mcroni'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = MCRONI_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    fa: Channel<Tuple<Map, Path>> = MCRONI_MODULE.out.fa

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        MCRONI_MODULE.out.fa
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        MCRONI_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
