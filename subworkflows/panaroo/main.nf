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
 * @modules panaroo_run
 *
 * @input gff
 * Channel containing gff data
 *
 * @output csv          Csv
 * @output aln          Aln
 * @output filtered_aln Filtered Aln
 * @output panaroo_csv  Panaroo Csv
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution logs from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PANAROO_RUN  } from '../../modules/panaroo/run/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow PANAROO {
    take:
    gff: Channel<Tuple<Map, Set<Path>>>

    main:
    PANAROO_RUN(gff)

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.csv
    aln: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.aln
    filtered_aln: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.filtered_aln
    panaroo_csv: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.panaroo_csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PANAROO_RUN.out.supplemental,
        PANAROO_RUN.out.csv,
        PANAROO_RUN.out.aln,
        PANAROO_RUN.out.filtered_aln,
        PANAROO_RUN.out.panaroo_csv
    ])
    logs: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.logs
    nf_logs: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.nf_logs
    versions: Channel<Tuple<Map, Path>> = PANAROO_RUN.out.versions
}
