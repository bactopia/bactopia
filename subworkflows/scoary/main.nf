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
 * @modules scoary as scoary_module
 *
 * @input csv
 * Channel containing csv data
 *
 * @input traits
 * Channel containing traits data
 *
 * @output csv      Csv
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SCOARY as SCOARY_MODULE } from '../../modules/scoary/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SCOARY {
    take:
    csv: Channel<Tuple<Map, Set<Path>>>
    traits: Path?

    main:    
    SCOARY_MODULE(csv, traits)

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.csv
    logs: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = SCOARY_MODULE.out.versions
}
