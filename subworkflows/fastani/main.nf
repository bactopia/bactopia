/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules csvtk_concat, fastani as fastani_module
 *
 * @input query
 * Channel containing query data
 *
 * @input reference
 * Channel containing reference data
 *
 * @output tsv        Tsv
 * @output merged_tsv Merged Tsv
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution logs from all processes
 * @output versions   Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { FASTANI as FASTANI_MODULE } from '../../modules/fastani/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow FASTANI {
    take:
    query: Channel<Tuple<Map, Set<Path>>>
    reference: Channel<Tuple<Map, Set<Path>>>

    main:
    FASTANI_MODULE(gather(query, 'query', 'fasta'), gather(reference, 'reference', 'fasta'))
    CSVTK_CONCAT(gather(FASTANI_MODULE.out.tsv, 'fastani'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = FASTANI_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
