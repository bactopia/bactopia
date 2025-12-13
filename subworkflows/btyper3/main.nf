/**
 * In silico taxonomic classification of Bacillus cereus group genomes.
 *
 * This subworkflow performs taxonomic classification of Bacillus cereus group
 * genomes using [BTyper3](https://github.com/lmc297/BTyper3), which provides
 * comprehensive classification including species, lineage, and toxin gene detection.
 * The results from individual samples are aggregated into a combined summary file.
 *
 * @status stable
 * @keywords Bacillus, cereus, taxonomy, typing, toxin genes
 * @tags complexity:simple input-type:single output-type:multiple features:aggregation
 * @citation btyper3
 *
 * @modules csvtk_concat, btyper3
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembly files in FASTA format for Bacillus cereus group classification
 *
 * @output tsv         BTyper3 classification results with detailed taxonomic information
 * @output merged_tsv  Combined TSV file containing classification results from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { BTYPER3 as BTYPER3_MODULE } from '../../modules/btyper3/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow BTYPER3 {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    BTYPER3_MODULE(assembly)
    CSVTK_CONCAT(gather(BTYPER3_MODULE.out.tsv, 'btyper3'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = BTYPER3_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.tsv,
        BTYPER3_MODULE.out.supplemental,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        BTYPER3_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
