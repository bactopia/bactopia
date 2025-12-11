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
 * @modules pirate as pirate_module
 *
 * @input gff
 * Channel containing gff data
 *
 * @output aln      Aln
 * @output csv      Csv
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PIRATE as PIRATE_MODULE } from '../../modules/pirate/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PIRATE {
    take:
    gff: Channel<Tuple<Map, Set<Path>>>

    main:
    PIRATE_MODULE(gff)

    emit:
    // Individual outputs
    aln: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.aln
    csv: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PIRATE_MODULE.out.supplemental,
        PIRATE_MODULE.out.aln,
        PIRATE_MODULE.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = PIRATE_MODULE.out.versions
}
