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
 * @modules snpdists as snpdists_module
 *
 * @input alignment
 * Channel containing alignment data
 *
 * @output tsv      Tsv
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SNPDISTS as SNPDISTS_MODULE } from '../../modules/snpdists/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SNPDISTS {
    take:
    alignment: Channel<Tuple<Map, Set<Path>>>

    main:
    SNPDISTS_MODULE(alignment)

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.tsv
    logs: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SNPDISTS_MODULE.out.versions
}
